#! /usr/bin/perl -w

use strict;
use Getopt::Long;

use strict;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use File::Path qw(make_path remove_tree);
use File::Basename;
use PerlIO::gzip;
use Compress::Zlib;
use Storable;
use IPC::Open2;
use FileHandle;

use constant FALSE => 0;
use constant TRUE  => 1;

my ($tarFile, $datFile, $desIdx, $splitChar, $maxCnt, $printHeader,
    $FR,$TB,$CO);

# $tarFile = "/Users/bbarnes/Documents/Projects/methylation/UCLA/workspace/orders/UCLA/build.9-21-2018/inputs/cpg-input.sorted.txt";
# $tarFile = "/Users/bbarnes/Documents/Projects/methylation/UCLA/workspace/orders/UCLA/build.9-23-2018/inputs/cpg-input.sorted.txt";

# $tarFile = "/Users/bbarnes/Documents/Projects/methylation/UCLA/input.match.bed";
# $datFile = "/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput/hs_ref_both_both_RefSeq37_probeDesignOutput.tsv.gz";

# $tarFile = '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/Redesign/data/LowIntProbes_Replacements.cgn-sorted.txt';
# $datFile = '/Users/bbarnes/Documents/Projects/workhorse/designs/mm10.improbeDesignOutput.fwdSeq-probes-scores-strands.tsv.gz';

$desIdx  = 0;
$printHeader = 0;

GetOptions("t=s"  => \$tarFile,
	   "d=s"  => \$datFile,

	   "idx=i" => \$desIdx,
	   "spl=s" => \$splitChar,

	   "FR=s"  => \$FR,
	   "TB=s"  => \$TB,
	   "CO=s"  => \$CO,

	   "header" => \$printHeader,

	   "max=i" => \$maxCnt,
    );

if (! defined($tarFile) || ! defined($datFile) || ! defined($desIdx)) {
    print STDERR "[Usage]: $0 -t tarFile -d datFile -idx desIdx [options]\n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n\n";
    exit(1);
}

my ($tarFH);
my ($tarFileName, $tarFilePath, $tarFileSuffix, $isTarBin)=openFile(\$tarFH,$tarFile);

my %tar = ();
while(<$tarFH>) {
    s/^\s+//;
    s/\s+$//;
    next if /^(\s)*$/;
    next if /^Seq_ID/;
    next if /^Probe_ID/;
    # next if /^Chromosome/;

    my $line = $_;
    my @data = split(/\t/, $line);

    my $id = $data[0];
    $tar{$id} = $id;

    # print STDERR "Target=$id\n";
}
close($tarFH);

print STDERR "[Targets]: Loaded " . scalar(keys %tar) . "\n\n";
print STDERR "[Loading]: Max Count=$maxCnt\n" if (defined($maxCnt));
print STDERR "[Loading]: dataFile=$datFile\n";

my ($datFH);
my ($datFileName, $datFilePath, $datFileSuffix, $isDatBin)=openFile(\$datFH,$datFile);

my $matCnt = 0;
my $misCnt = 0;
my $totCnt = 0;
while(<$datFH>) {
    s/^\s+//;
    s/\s+$//;
    # next if /^Seq_ID/;
    
    my $line = $_;
    if ($printHeader && $line =~ m/^Seq_ID/) {
	print $line."\n" if ($printHeader);
	next;
    }

    my @data = split(/\t/, $line);

    # print STDERR "Line=$line\n";

    # print STDERR "LINE=$line\n";
    if (!defined($line) || length($line)==0) {
	die("[ERROR]: line is null!\n");
    }

    my $cpg = $data[$desIdx];
    my $len = length($cpg);

    # if (!defined($len) || $len != 24) {
    if (!defined($len)) {
	# next;
	print STDERR "\n\n";
	print STDERR "TotCnt=$totCnt\n";
	print STDERR "Line='$line'\n";
	print STDERR "DataCpg=$cpg, len=$len\n";
	print STDERR "Data=@data\n";
	exit(1);
    }

    my @cpgDat;
    if (defined($splitChar)) {
	@cpgDat = split($splitChar, $cpg);
	$cpg=$cpgDat[0];
    }

    my $prbFR=$data[20];
    my $prbTB=substr($data[21],0,1);
    my $prbCO=$data[22];

    # print STDERR "$cpg\tFR=$prbFR, TB=$prbTB, CO=$prbCO\n";

    # next if (defined($FR) && defined($data[20]) && $data[21] ne $FR);
    # next if (defined($TB) && defined($data[21]) && $data[22] ne substr($TB, 0,1));
    # next if (defined($CO) && defined($data[22]) && $data[23] ne $CO);

    # $cpg=$cpg."_".$prbTB.$prbCO;

    if (defined($tar{$cpg})) {
	$tar{$cpg}++;
	print "$line\n";
	$matCnt++;
    } else {
	$misCnt++;
    }

    print STDERR "[Status]: mat=$matCnt, mis=$misCnt, tot=$totCnt; desIdx=$desIdx; cpg=$cpg\n" if ($totCnt % 100000 == 0);
    # print STDERR "[Status]: mat=$matCnt, mis=$misCnt, tot=$totCnt; desIdx=$desIdx; cpg=$cpg\n\tDAT='@data'\n" if ($totCnt % 1000000 == 0);

    last if (defined($maxCnt) && $totCnt >= $maxCnt);
    $totCnt++;
}
close($datFH);
print STDERR "[Done]: mat=$matCnt, mis=$misCnt, tot=$totCnt\n";



#
# Subroutines::
#

sub trimSuffix {
    my $str=shift;

    for (my $ii = 0; $ii < 5; $ii++) {
	$$str=~ s/\.gz$//;  $$str=~ s/\.fa$//;
        $$str=~ s/\.bed$//; $$str=~ s/\.txt$//; $$str=~ s/\.csv$//; $$str=~ s/\.tsv$//;
        $$str=~ s/\.bam$//; $$str=~ s/\.sam$//; $$str=~ s/\.vcf$//; $$str=~ s/\.fas$//;
        $$str=~ s/\.bsp$//; $$str=~ s/\.manifest$//; $$str=~ s/\.manifests$//;
        $$str=~ s/\.sorted$//; $$str =~ s/\.merged$//; $$str =~ s/\.intersect$//;
	$$str=~ s/\.[^\.]-sorted$//;
    }
}

sub openFile {
    my $func=(caller(0))[3];
    my $fh=shift;
    my $file=shift;

    die("[$func]: ERROR File does not exist '$file'\n")
        if (! defined($file) || ! -e $file);

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    my $isBin = ($fileName =~ /\.bin$/) ? 1 : 0;
    trimSuffix(\$fileName);

    $$fh = gzopen($file, "rb") or die("Failed to open binary gzip stream=$file for reading: $!")
        if ($isBin && $fileSuffix eq ".gz");

    open($$fh, "<:gzip", $file) or die("Failed to open gzip stream=$file for reading: $!")
        if (! $isBin && $fileSuffix eq ".gz");

    open($$fh, "<:raw", $file) or die("Failed to open binary file=$file for reading: $!")
        if ($isBin && $fileSuffix ne ".gz");

    open($$fh, "<$file") or die("Failed to open regular file=$file for reading: $!")
	if (! $isBin && $fileSuffix ne ".gz");

    return(($fileName, $filePath, $fileSuffix, $isBin));
}


## End of file
