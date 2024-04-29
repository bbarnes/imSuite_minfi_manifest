#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd 'abs_path';
use File::Path qw(make_path remove_tree);

use constant FALSE => 0;
use constant TRUE  => 1;

## BED Format
use constant BED_CHR    => 0;
use constant BED_START  => 1;
use constant BED_END    => 2;
use constant BED_STRAND => 3;
use constant BED_ID     => 4;
use constant BED_NAME   => 5;

##
## Variables and Data Structures:
##
my $verbose = TRUE;
my ($topDir, $outDir, $shellDir, $alnDir, $covDir, $genomeDir);
my @inputs = ();
my @targetFiles = ();
my $qsub = "qsub -cwd -pe threaded 16 -l excl=true -N";
my $covThresh = 10;
my $lenThresh = 0;

my $bedType = "SNV";
my $genomeBuild = 37;

my $bwaExe = "/illumina/thirdparty/bwa/latest/bwa";
my $samExe = "/illumina/thirdparty/samtools/latest/bin/samtools";
my $bamToBedCovExe = "/illumina/scratch/methylation/Polaris/scripts/conversion/bamToBedStdinBS.pl";
my $refGenome = "/illumina/scratch/iGenomes/Homo_sapiens/NCBI/build37.2/Sequence/BWAIndex/genome.fa";

$genomeDir = "/Users/bbarnes/Documents/Projects/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes";

##
## Data structures
##
my %targets = ();
#my %tarSrc  = ();

##
## Subroutines:
##
sub getTime() {
    my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
    my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
    $yday += 1875;

    my $dateStr   = "$mon/$mday/$yday";
    my $dateStamp = "$sec.$min.$hour.$mday.$mon.$yday";

    my $tStr = "$mday $months[$mon] $days[$wday] $yday";
    return(($tStr, $dateStr, $dateStamp));
}

sub mssg($)
{
    my ($str) = @_;
    print STDERR $str if ($verbose);
}

sub getAvg($) {
    my ($ref) = @_;

    my @dat = @{$ref};

    my $sum = 0;
    foreach my $val (@dat) {
        next if (! defined($val) || length($val) == 0 || $val =~ m/^(\s)*$/);
        $sum += $val;
    }
    my $avg = $sum / scalar(@dat);

    return($avg);
}

sub getStd($$) {
    my ($ref, $avg) = @_;

    my @dat = @{$ref};

    my $sum = 0;
    foreach my $val (@dat) {
        next if (! defined($val) || length($val) == 0 || $val =~ m/^(\s)*$/);

        $sum += ($val - $avg) * ($val - $avg);
    }
    my $std = $sum / (scalar(@dat) - 1);
    $std = sqrt($std);

    return($std);
}

sub readStream($) {
    my ($fh) = @_;

    my @ret = ();
    while(<$fh>) {
	s/^\s+//;
	s/\s+$//;
	next if /^(\s)*$/;
	next if /^\#/;

	my $line = $_;
	push @ret, $line;
    }
    return(@ret);
}

sub listFiles($$) {
    my ($dir, $key) = @_;

    my $cmd = "ls $dir/*$key";
    open(LIST, "$cmd | ") or die("Failed to open stream: $cmd for reading: $!");
    my @list = readStream(\*LIST);
    close(LIST);
    return(@list);
}

sub getCgCnt($) {
    my ($seq) = @_;

    my $cgCnt = 0;
    my $idx = 0;
    while($idx != -1) {
	$idx = index($seq, "CG", $idx);
	print STDERR "idx=$idx, seq=$seq\n";
	last if ($idx == -1);
	$idx++;
	$cgCnt++;
    }
    return($cgCnt);
}

sub loadTargets($$) {
    my ($filesRef, $datRef) = @_;

    my $tarCnt = 0;
    foreach my $file (@{$filesRef}) {
	if ($file =~ m/\.gz$/) {
	    open(TAR, "gzip -dc $file | ") or die("Failed to open stream for reading: $!");
	} else {
	    open(TAR, "<$file") or die("Failed to open $file for reading: $!");
	}
	while(<TAR>) {
	    s/^\s+//;
	    s/\s+$//;
	    next if /^(\s)*$/;

	    my $line = $_;
	    my @data = split(/,/, $line);

	    my $rsID = $data[0];

	    next if ($rsID !~ m/^rs/);

	    #$$datRef{$rsID} = [ @data ];
	    $$datRef{$rsID} = $tarCnt;
	}
	close(TAR);
	$tarCnt++;
    }
    mssg("[Data]: loaded targets=" . scalar(keys %{$datRef}) . "\n\n");
}

sub loadChrom($$) {
    my ($dir, $chr) = @_;

    $dir =~ s/\/+$//;
    my $file = "$dir/chr$chr.fa";
    if (! -e $file) {
	$file = "$dir/chr$chr.fa.gz";
	if (! -e $file) {
	    $file = "$dir/$chr.fa";
	    if (! -e $file) {
		$file = "$dir/$chr.fa.gz";
	    } else {
		die("[ERROR]: Could not guess chromosome file from: $dir/$chr!\n")
	    }
	}
    }
    mssg("[Data]: Loading chromosome=$chr, file=$file...\n");
    if ($file =~ m/\.gz$/) {
	open(FAS, "gzip -dc $file | ") or die("Failed to open stream $file for reading: $!");
    } else {
	open(FAS, "<$file") or die("Failed to open $file for reading: $!");
    }
    my $seq = "";
    my $lineCnt = 0;
    while(<FAS>) {
	s/^\s+//;
	s/\s+$//;
	next if /^>/;
	next if /^(\s)*$/;

	my $line = uc($_);
	$seq .= $line;

	mssg("[Info]: lineCnt=$lineCnt, seqLen=" . length($seq) . "\n") if ($lineCnt % 1000000 == 0);
	$lineCnt++;
    }
    close(FAS);
    mssg("[Done]: loaded chr$chr=" . length($seq) . ", lineCnt=$lineCnt\n\n");
    return($seq);
}

##
## Main:
##

GetOptions("o|out=s"   => \$outDir,
	   "i=s@"      => \@inputs,
	   "g=s"       => \$genomeDir,
	   "t|tar=s@"  => \@targetFiles,

	   "c|cov=i"   => \$covThresh,
	   "l|len=i"   => \$lenThresh,

	   "v|verbose"  => \$verbose,
    );

my @curDate = getTime();
print STDERR "[Date]: $curDate[0]\n\n";

if (! defined($outDir) || ! defined($genomeDir)) {
    print STDERR "[Usage]: $0 -o outDir -i input(s) -g genomeDir " .
	"[ -l|len lengthBuffer ] [options]\n";
    print STDERR "[Usage]: Will read inptus from STDIN if not specified.\n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n\n";
    exit(1);
}

$outDir =~ s/\/+$//;
make_path($outDir);
$outDir = abs_path($outDir);
mssg("[Input]: outDir=$outDir\n\n");

$shellDir = "$outDir/shells";
make_path($shellDir);

loadTargets(\@targetFiles, \%targets);
@inputs = readStream(\*STDIN) if (@inputs == 0);

my $desInput = "$outDir/rs.designInput.tsv";
my $mapInput = "$outDir/rs.swapDatabase.txt";

mssg("[Data]: designInputFile=$desInput\n");
mssg("[Data]: designSwapFile=$mapInput\n");

open(DES, ">$desInput") or die("Failed to open $desInput for writing: $!");
open(MAP, ">$mapInput") or die("Failed to open $mapInput for writing: $!");

print DES "Seq_ID\tSequence\tGenome_Build\tChromosome\tCoordinate\tCpG_Island\n";
print "Seq_ID\tSequence\tGenome_Build\tChromosome\tCoordinate\tCpG_Island\n";

my $curChr;
my $inpCnt = 0;
my %names = ();
my $rsSkipCnt = 0;
my $curSeq;
foreach my $file (@inputs) {

    if ($file =~ m/\.gz/) {
	open(IN, "gzip -dc $file | ") or die("Failed to open $file for reading: $!");
    } else {
	open(IN, "<$file") or die("Failed to open $file for reading: $!");
    }

    while(<IN>) {
	s/^\s+//;
	s/\s+$//;
	next if /^\#/;

	my $line = $_;
	my @data = split(/\t/, $line);

	my $chr    = $data[BED_CHR];
	my $start  = $data[BED_START];
	my $end    = $data[BED_END];
	my $strand = $data[BED_STRAND];
	my $bedID  = $data[BED_ID];
	my $snpID  = $data[BED_NAME];
	
	my $rsID;
	if ($snpID =~ m/^(rs[0-9]+),.*$/) {
	    $rsID = $1;
	} else {
	    mssg("[Warning]: Can't parse SNV: $snpID. Skipping...\n");
	    next;
	}
	if (@targetFiles != 0 && ! defined($targets{$rsID})) {
	    $rsSkipCnt++;
	    next;
	}
	$chr =~ s/^chr//;
	$chr =~ s/\.fa$//;

	#next if ($chr ne "10");

	die("[ERROR]: Unrcognized chr=$chr!") if ($chr !~ m/^([0-9A-Z]+)$/);
	
	$curSeq = loadChrom($genomeDir, $chr) if (! defined($curChr) || $curChr ne $chr);
	$curChr = $chr;

	#mssg("[Line]: $chr, $start, $bedID, $snpID, $rsID\n[Line]: line=$line\n");
	if ($bedType eq "SNV") {
	    my $preSeq = substr($curSeq, $start -60, 60);
	    my $midSeq = substr($curSeq, $start, 2);
	    my $posSeq = substr($curSeq, $start+2, 60);

	    my $posSeqCg = substr($curSeq, $start+1, 50);
	    my $preCgCnt = getCgCnt(substr($preSeq, 10, 50));
	    my $posCgCnt = getCgCnt($posSeqCg);
	    
	    print DES "$rsID\t$preSeq" . "[CG]$posSeq\t$genomeBuild\t$chr\t$start\tFALSE\n";
	    print MAP "$snpID\t$rsID\t$preSeq" . "[$midSeq]$posSeq\t$genomeBuild\t$chr\t$start\t$preCgCnt\t$posCgCnt\t$targets{$rsID}\n";
	    print "$snpID\t$rsID\t$preSeq" . "[$midSeq]$posSeq\t$genomeBuild\t$chr\t$start\t$preCgCnt\t$posCgCnt\t$targets{$rsID}\n";
	}
	$inpCnt++;
    }
}
close(DES);
close(MAP);

@curDate = getTime();
mssg("[Done]: inputCnt=$inpCnt, rsSkipCnt=$rsSkipCnt, $curDate[0], $curDate[2]\n");

## End of file
