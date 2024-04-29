#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd;
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
$topDir = abs_path(getcwd);
my $jobName = "job";

## Input Data Structures:
my @inpFiles = ();

my @srcNames = ();
my @dmrFiles = ();

my @bamFiles = ();
my @desFiles = ();
my @annFiles = ();

my %srcNam  = ();
my %dmrDat  = ();
my %bamDat  = ();
my %desDat  = ();
my %annDat  = ();

##
## Program Parameters
##
my ($dmrMax, $selMax, $selSweep);
my $statusCnt = 10000;
my $inputMax;

##
## Cluster/Mac Standard Files:
##
my $qsub = "qsub -cwd -pe threaded 16 -l excl=true -N";
my ($bwaExe, $samExe, $bamToBedCovExe, $refGenome, $hg19nc, $hg19, $chrDir);
if ($topDir =~ m/^\/illumina\//) {
    $bwaExe = "/illumina/thirdparty/bwa/latest/bwa";
    $samExe = "/illumina/thirdparty/samtools/latest/bin/samtools";
    $bamToBedCovExe = "/illumina/scratch/methylation/Polaris/scripts/conversion/bamToBedStdinBS.pl";
    $refGenome = "/illumina/scratch/iGenomes/Homo_sapiens/NCBI/build37.2/Sequence/BWAIndex/genome.fa";
} else {
    $hg19   = "/Users/bbarnes/Documents/Projects/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/noChr.genome.fa";
    $hg19nc = "/Users/bbarnes/Documents/Projects/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa";
    $chrDir = "/Users/bbarnes/Documents/Projects/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes";
}

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

sub trimSuffix($) {
    my ($fileName) = @_;

    $fileName =~ s/\.gz$//;
    $fileName =~ s/\.bed$//;
    $fileName =~ s/\.sorted$//;
    $fileName =~ s/\.intersect$//;
    $fileName =~ s/\.merged$//;
    $fileName =~ s/\.log$//;
    $fileName =~ s/\.csv$//;
    $fileName =~ s/\.txt$//;
    $fileName =~ s/\.tsv$//;
    $fileName =~ s/\.bam$//;
    $fileName =~ s/\.sam$//;
    $fileName =~ s/\.vcf$//;
    $fileName =~ s/\.gz$//;
    $fileName =~ s/\.bed$//;
    $fileName =~ s/\.sorted$//;
    $fileName =~ s/\.intersect$//;
    $fileName =~ s/\.merged$//;
    $fileName =~ s/\.log$//;
    $fileName =~ s/\.csv$//;
    $fileName =~ s/\.txt$//;
    $fileName =~ s/\.tsv$//;
    $fileName =~ s/\.bam$//;
    $fileName =~ s/\.sam$//;
    $fileName =~ s/\.vcf$//;

    return($fileName);
}

##
## Main:
##

GetOptions("o|out=s"     => \$outDir,

	   "src=s@  "    => \@srcNames,
	   "dmr=s@"      => \@dmrFiles,
	   "des=s@"      => \@desFiles,
	   "b|bamss@"    => \@bamFiles,
	   "a|ann=s@"    => \@annFiles,

	   "dmrMax=i"    => \$dmrMax,
	   "selMax=i"    => \$selMax,

	   "sweep"       => \$selSweep,

	   "v|verbose"  => \$verbose,
    );

my @curDate = getTime();
print STDERR "[Date]: $curDate[0]\n\n";

if (! defined($outDir)) {
    print STDERR "[Usage]: $0 -o outDir "
	. " [ -src sourcName(s) ]\n"
	. " [ -dmr dmrBed(s) ]\n"
	. " [ -des desBed(s) ]\n"
	. " [ -b|bam bam(s) ]\n"
	. " [ -a|ann annotationBed(s) ]\n"
	

	. " [ options ]\n";
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

@inpFiles = readStream(\*STDIN) if (@inpFiles == 0);

my $inpCnt = 0;
my %bamNames = ();
foreach my $file (@inpFiles) {
    $fileFile = abs_path($fileFile);

    my $fileName;
    if ($fileFile =~ m/^.*\/([^\/]+)$/) {
	$fileName = $1;
	$fileName = trimSuffix($fileName);
    }
    if (defined($fileNames{$fileName})) {
	mssg("[Warning]: skipping name=$fileName, multiple definitions!\n");
	next;
    }
    $fileNames{$fileName} = 1;

    my %bedNames = ();
    foreach my $bedFile (@bedFiles) {
	my $bedName;
	if ($bedFile =~ m/^.*\/([^\/]+)$/) {
	    $bedName = $1;
	    $bedName = trimSuffix($bedName);
	}
	if (defined($bedNames{$bedName})) {
	    mssg("[Warning]: skipping name=$bedName, multiple definitions!\n");
	    next;
	}
	$bedNames{$bedName} = 1;

	my $intFile = "$outDir/$fileName.$bedName.bedCov.vcf";
	$cmd = "samtools mpileup -u -v --output-MQ -t ADF,ADR -f $hg19 -l $bedFile $fileFile > $intFile";
	my $retVal = system($cmd);
	mssg("[mpileup]: $cmd...\n");

	die("[ERROR]: mpileup command failed: '$cmd'") if ($retVal);
	last;
	if ($intFile=~ m/\.gz$/) {
	    open(IN, "gzip -dc $intFile | ")or die("Failed to open $intFile STREAM for reading: $!");
	} else {
	    open(IN, "<$intFile") or die("Failed to open $intFile for reading: $!");
	}

	my $misCnt  = 0;
	my $matCnt  = 0;
	my $totCnt  = 0;

	my $lineCnt = 0;
	while(<IN>) {
	    s/^\s+//;
	    s/\s+$//;

	    my $line = $_;
	    my @data = split(/\t/, $line);

	    my $chr   = $data[VCF_CHROM];
	    my $pos   = $data[VCF_POS];

	    $lineCnt++;
	    last if (defined($inputMax) && $lineCnt >= $inputMax);
	}
	close(IN);
	mssg("[Input]: misCnt=$misCnt, matCnt=$matCnt, totCnt=$totCnt\n");
	mssg("[Input]: lineCnt=$lineCnt\n\n");
	
	next;
	my $outFile   = "$bedDir/$bedName.counts.bed";
	my $sortedBed = "$bedDir/$bedName.sorted.bed";
	open(OUT, ">$outFile") or die("Failed to open $outFile for writing: $!");
	
	close(OUT);
	
	$cmd = "sort -n -k 2,2 $outFile | sort -s -k 1,1 > $sortedBed";
	system($cmd);
	system("rm -f $sortedBed.gz") if (-e "$sortedBed.gz");
	system("rm -f $sortedBed.gz.tbi") if (-e "$sortedBed.gz.tbi");
	system("bgzip $sortedBed");
	system("tabix $sortedBed.gz");

    next;
    my $runShell = "$shellDir/run.$bedName.$inpCnt.sh";

    open(SHELL, ">$runShell") or die("Failed to open $runShell for writing: $!");
    print SHELL "cd $outDir\n\n";
    print SHELL "touch $outDir/run.$bedName.$inpCnt.finished.txt";
    close(SHELL);

    my $launchShell = "$shellDir/launch.$bedName.$inpCnt.sh";
    open(SHELL, ">$launchShell") or die("Failed to open $launchShell for writing: $!");
    print SHELL "$runShell\n";
    close(SHELL);

    system("chmod 777 $runShell");
    system("chmod 777 $launchShell");

    $cmd = "$qsub $jobName$inpCnt $launchShell";
    mssg("[qsub]: '$cmd'\n");
    system($cmd);
    }
}

@curDate = getTime();
mssg("[Done]: inputCnt=$inpCnt, $curDate[0], $curDate[2]\n");

## End of file
