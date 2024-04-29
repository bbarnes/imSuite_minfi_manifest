#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use File::Path qw(make_path remove_tree);
use File::Basename;

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

## BAM Format:
use constant BAM_QNAME    => 0;
use constant BAM_FLAG     => 1;
use constant BAM_RNAME    => 2;
use constant BAM_POS      => 3;
use constant BAM_MAPQ     => 4;
use constant BAM_CIGAR    => 5;
use constant BAM_RNEXT    => 6;
use constant BAM_PNEXT    => 7;
use constant BAM_TLEN     => 8;
use constant BAM_SEQ      => 9;
use constant BAM_QUAL     => 10;
use constant BAM_SAMPLE   => 11;

##
## Variables and Data Structures:
##
my $verbose = TRUE;
my ($topDir, $outDir, $logDir, $bedDir);
$topDir = abs_path(getcwd);
my $jobName = "job";

##
## Program Parameters
##
my $statusCnt = 10000;
my $inpMax;

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

sub mssg
{
    my $str = shift;
    my $level = shift;

    if (defined($level)) {
	print STDERR $str if ($verbose >= $level);
    } else {
	print STDERR $str if ($verbose);
    }
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
    my ($str) = @_;

    for (my $ii = 0; $ii < 5; $ii++) {
    $str =~ s/\.gz//; $str =~ s/\.fa//;
    $str =~ s/\.bed//; $str =~ s/\.txt//; $str =~ s/\.csv//; $str =~ s/\.tsv//;
    $str =~ s/\.bam//; $str =~ s/\.sam//; $str =~ s/\.vcf//; $str =~ s/\.fas//;
    $str =~ s/\.sorted//; $str =~ s/\.merged//; $str =~ s/\.intersect//;
    }
    return($str);
}

sub bgzipIdx {
    my ($file, $clean) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    my $sortedBed = "$filePath$fileName.sorted.bed";

    my $cmd = "sort -k 1,1V -k 2,2n -k 3,3n $file > $sortedBed";
    die("[bgzipIdx]: Failed sort: '$cmd'") if (system("$cmd"));
    system("rm -f $sortedBed.gz") if (-e "$sortedBed.gz");
    system("rm -f $sortedBed.gz.tbi") if (-e "$sortedBed.gz.tbi");
    die("[BgzipIdx]: Failed bgzip $sortedBed") if (system("bgzip $sortedBed"));
    $sortedBed .= ".gz";
    die("[BgzipIdx]: Failed tabix $sortedBed") if (system("tabix $sortedBed"));
    system("rm -f $file") if ($clean);
    return($sortedBed);
}

sub samToBam {
    my ($file, $clean) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    my $bam = "$filePath$fileName.bam";
    my $sortedBam = "$filePath$fileName.sorted.bam";

    my $cmd = "samtools view -S -b $file > $bam";
    die("[samToBam]: Failed sam <-> bam: $cmd") if (system($cmd));
    $cmd = "samtools sort $bam > $sortedBam";
    die("[samToBam]: Failed sort: $cmd") if (system($cmd));
    $cmd = "samtools index $sortedBam";
    die("[samToBam]: Failed index: $cmd") if (system($cmd));
    system("rm -f $file") if ($clean);
    system("rm -f $bam") if ($clean);
}

##
## Program Subroutines:
##
sub loadBed {
    my ($file, $datRef, $bedMax, $sortRef, $writeBed) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    $fileName = trimSuffix($fileName);

    my ($sortedBed);
    my $outBed = "$bedDir/$fileName.bed";
    open(OUT, ">$outBed") or die("Failed to open $outBed for writing: $!") if ($writeBed);

    if ($fileSuffix eq ".gz") { open(BED, "gzip -dc $file | ") or die("Failed to open $file stream: $!") }
    else { open(BED, "<$file") or die("Failed to open $file: $!") }

    my %keys = ();
    my $inpCnt = 0;
    while(<BED>) {
	s/^\s+//;
	s/\s+$//;

	my $line = $_;
	my @data = split(/\t/, $line);

	my $datKey = $data[BED_NAME];

	## Check for Duplicates
	die("[ERROR]: Multiple key definitions for $datKey!") if (defined($keys{$datKey}));
	$keys{$datKey} = 1;

	## Store Data
	$$datRef{$datKey} = [ @data ];
	push @{$sortRef}, $datKey if (defined($sortRef));

	## Write Data
	print OUT "$line\n" if ($writeBed);
	    
	## Status
	mssg("[Bed]: statusCnt=$inpCnt...\n", 4)
	    if (defined($statusCnt) && $inpCnt % $statusCnt == 0);
	
	$inpCnt++;
	last if (defined($bedMax) && $inpCnt >= $bedMax);
    }
    close(BED);
    if ($writeBed) {
	close(BED);
	$sortedBed = bgzipIdx($outBed);
    }
    mssg("[Bed]: Loaded: " . scalar(keys %keys) . " targets\n", 3);

    return($sortedBed) if ($writeBed);
}


##
## Main:
##

GetOptions("o|out=s"     => \$outDir,

	   "v|verbose=i" => \$verbose,
    );

my @curDate = getTime();
print STDERR "[Date]: $curDate[0]\n\n";

if (! defined($outDir)) {
    print STDERR "[Usage]: $0 -o outDir "
	. " [ v|verbose verbosityLevel ]\n"
	. " [ options ]\n\n";
    print STDERR "[Usage]: Will read inptus from STDIN if not specified.\n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n\n";
    exit(1);
}

$outDir =~ s/\/+$//;
make_path($outDir);
$outDir = abs_path($outDir);
mssg("[Input]: outDir=$outDir\n\n");

$logDir = "$outDir/log";
make_path($logDir) if (! -e $logDir);
system("rm -f $logDir/* &> $outDir/cleanLog.log");

$bedDir = "$outDir/bed";
make_path($bedDir) if (! -e $bedDir);
system("rm -f $bedDir/* &> $bedDir/cleanBed.log");

##
## Program:
##

my $inpCnt=0;


##
## Done!
##
@curDate = getTime();
mssg("[Done]: inputCnt=$inpCnt, $curDate[0], $curDate[2]\n\n");

## End of file
