#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd 'abs_path';
use File::Path qw(make_path remove_tree);

use constant FALSE => 0;
use constant TRUE  => 1;

## BED Format:
use constant BED_CHR    => 0;
use constant BED_START  => 1;
use constant BED_END    => 2;
use constant BED_NAME   => 3;
use constant BED_SCORE  => 4;
use constant BED_STRAND => 5;

use constant BED_OFFSET => 5;

## VCF Format:
use constant VCF_CHROM   => 0;
use constant VCF_POS     => 1;
use constant VCF_ID      => 2;
use constant VCF_REF     => 3;
use constant VCF_ALT     => 4;
use constant VCF_QUAL    => 5;
use constant VCF_FILTER  => 6;
use constant VCF_INFO    => 7;
use constant VCF_FORMAT  => 8;
use constant VCF_SAMPLE  => 9;
use constant VCF_GT      => 10;
use constant VCF_QC      => 11;
use constant VCF_REF_CNT => 12;
use constant VCF_ALT_CNT => 13;

##
## Variables and Data Structures:
##
my $verbose = TRUE;
my ($topDir, $outDir, $shellDir, $logDir, $bedDir, $datDir, $cmd);
my @inputs = ();
my $jobName = "job";

my $statusCnt = 10000;
my $inputMax;

my $refOffset  = BED_OFFSET;

my $qsub = "qsub -cwd -pe threaded 16 -l excl=true -N";
my $bwaExe = "/illumina/thirdparty/bwa/latest/bwa";
my $samExe = "/illumina/thirdparty/samtools/latest/bin/samtools";
my $bamToBedCovExe = "/illumina/scratch/methylation/Polaris/scripts/conversion/bamToBedStdinBS.pl";
my $refGenome = "/illumina/scratch/iGenomes/Homo_sapiens/NCBI/build37.2/Sequence/BWAIndex/genome.fa";
my $hg19 = "/Users/bbarnes/Documents/Projects/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/noChr.genome.fa";

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

##
## Main:
##

GetOptions("o|out=s"   => \$outDir,
	   "i=s@"      => \@inputs,

	   "ro|refOff=i" => \$refOffset,

	   "v|verbose"  => \$verbose,
    );

my @curDate = getTime();
print STDERR "[Date]: $curDate[0]\n\n";

if (! defined($outDir)) {
    print STDERR "[Usage]: $0 -o outDir -i input(s) [ -ro|refOff refDataOffset ] "
	. " [options]\n";
    print STDERR "[Usage]: Will read inptus from STDIN if not specified.\n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n\n";
    exit(1);
}

$outDir =~ s/\/+$//;
make_path($outDir);
$outDir = abs_path($outDir);
mssg("[Input]: outDir=$outDir\n\n");

$shellDir = "$outDir/shells";
make_path($shellDir) if (! -e $shellDir);

$logDir = "$outDir/log";
make_path($logDir) if (! -e $logDir);

$bedDir = "$outDir/bed";
make_path($bedDir) if (! -e $bedDir);

$datDir = "$outDir/dat";
make_path($datDir) if (! -e $datDir);

@inputs = readStream(\*STDIN) if (@inputs == 0);

my $inpCnt = 0;
my %names = ();
foreach my $inFile (@inputs) {
    $inFile = abs_path($inFile);

    my $name = $inFile;
    if ($inFile =~ m/^.*\/([^\/]+)$/) {
	$name = $1;
    }
    if (defined($names{$name})) {
	mssg("[Warning]: skipping name=$name, multiple definitions!\n");
	next;
    }
    $names{$name} = 1;

    if ($inFile =~ m/\.gz$/) {
	open(IN, "gzip -dc $inFile | ") or die("Failed to open $inFile STREAM for reading: $!");
    } else {
	open(IN, "<$inFile") or die("Failed to open $inFile for reading: $!");
    }

    my $lineCnt = 0;
    my %bedDat = ();
    my %intDat = ();

    my $matCnt = 0;
    my $misCnt = 0;
    my $totCnt = 0;
    while(<IN>) {
	s/^\s+//;
	s/\s+$//;

	my $line = $_;
	my @data = split(/\t/, $line);

	my $refName  = $data[BED_NAME];
	my $canName  = $data[$refOffset+BED_NAME+1];
	my $canStart = $data[$refOffset+BED_START+1];
	my $canEnd   = $data[$refOffset+BED_END+1];

	$bedDat{$refName} = [ @data[BED_CHR...$refOffset] ] if (! defined($bedDat{$refName}));

	#mssg("$canStart\t$canEnd\n");
	if ($canStart eq "-1" && $canEnd eq "-1") {
	    $misCnt++;
	    next;
	}
	$intDat{$refName}{$canName} = 0 if (! defined($intDat{$refName}{$canName}));
	$intDat{$refName}{$canName}++;

	$lineCnt++;
        last if (defined($inputMax) && $lineCnt >= $inputMax);
    }
    close(IN);
    mssg("[Input]: misCnt=$misCnt, lineCnt=$lineCnt\n\n");

    ##
    ## Write new BED
    ##
    my $outFile   = "$bedDir/$name.counts.bed";
    my $sortedBed = "$bedDir/$name.sorted.bed";
    open(OUT, ">$outFile") or die("Failed to open $outFile for writing: $!");
    foreach my $refKey (sort keys %bedDat) {
	my $intCnt = 0;
	my $intStr = "";

	if (defined($intDat{$refKey})) {
	    $intCnt = scalar(keys %{$intDat{$refKey}});
	    foreach my $canKey (sort keys %{$intDat{$refKey}}) {
		$intStr .= "$canKey;";
	    }
	    $intStr =~ s/\;+$//;
	}

	print OUT join("\t", @{$bedDat{$refKey}}) . "\t$intCnt\t$intStr\n";
    }
    close(OUT);

    $cmd = "sort -n -k 2,2 $outFile | sort -s -k 1,1 > $sortedBed";
    system($cmd);
    system("rm -f $sortedBed.gz") if (-e "$sortedBed.gz");
    system("rm -f $sortedBed.gz.tbi") if (-e "$sortedBed.gz.tbi");
    system("bgzip $sortedBed");
    system("tabix $sortedBed.gz");

    next;
    my $runShell = "$shellDir/run.$name.$inpCnt.sh";

    open(SHELL, ">$runShell") or die("Failed to open $runShell for writing: $!");
    print SHELL "cd $outDir\n\n";
    print SHELL "touch $outDir/run.$name.$inpCnt.finished.txt";
    close(SHELL);

    my $launchShell = "$shellDir/launch.$name.$inpCnt.sh";
    open(SHELL, ">$launchShell") or die("Failed to open $launchShell for writing: $!");
    print SHELL "$runShell\n";
    close(SHELL);

    system("chmod 777 $runShell");
    system("chmod 777 $launchShell");

    $cmd = "$qsub $jobName$inpCnt $launchShell";
    mssg("[qsub]: '$cmd'\n");
    system($cmd);
}

@curDate = getTime();
mssg("[Done]: inputCnt=$inpCnt, $curDate[0], $curDate[2]\n");

## End of file
