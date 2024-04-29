#! /usr/bin/perl -w

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

my $srcDir;
my $outDir;
my $type="idats";

my %dat=();

GetOptions("s=s" => \$srcDir,
	   "o=s" => \$outDir,
	   "t=s" => \$type,
    );

if (! defined($outDir)) {
    print STDERR "[Usage]: $0 -t type -o outDir < find-list [options]\n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n";
    exit(1);
}

# die("[ERROR]: Source Directory must exist='$srcDir'\n\n") if (! -e $srcDir);
# $srcDir =~ s/\/+$//;
$outDir =~ s/\/+$//;
if (! -e $outDir) { make_path($outDir); }

my $totCnt = 0;
my $matCnt = 0;
my $misCnt = 0;
my $newCnt = 0;

while(<STDIN>) {
    s/^\s+//;
    s/\s+$//;

    my $file=$_;
    # $file =~ s/^\.\///;
    $file = abs_path($file);

    # $file = $srcDir."/".$file;
    die("[ERROR]: file=$file does not exist!\n") if (! -e $file);
    $totCnt++;

    my ($basecode, $poscode, $fileName, $color);

    my $isGZ = 0;
    $isGZ=1 if ($file =~ m/\.gz$/);

    if ($type eq "idats") {
	if ($file =~ m/^.*\/(([A-Z0-9]+)_([RC0-9]+)_(Grn|Red)\.idat$)/ ||
	    $file =~ m/^.*\/(([A-Z0-9]+)_([RC0-9]+)_(Grn|Red)\.idat.gz$)/ ||
	    $file =~ m/^.*\/.*_(([A-Z0-9]+)_([RC0-9]+)_(Grn|Red)\.idat$)/ ||
	    $file =~ m/^.*\/.*_(([A-Z0-9]+)_([RC0-9]+)_(Grn|Red)\.idat.gz$)/) {
	    $fileName=$1;
	    $basecode=$2;
	    $poscode=$3;
	    $color=$4;
	} else {
	    print STDERR "[Skipping]: Unrecognized IDAT file=$file\n";
	    next;
	}
    } else {
	if ($file =~ m/^.*\/(([A-Z0-9]+)_([RC0-9]+)_([0-9]+)\.dmap.gz)/) {
	    $fileName=$1;
	    $basecode=$2;
	    $poscode=$3;
	    $color=$4;
	} else {
	    print STDERR "[Skipping]: Unrecognized DMAP file=$file\n";
	    next;
	}
    }
    if (defined($dat{$fileName})) {
	print STDERR "[Duplicate]: $fileName\n";
	next;
    }

    # print STDERR "fileName=$file, basecode=$basecode, poscode=$poscode\n";

    my $distDir = $outDir."/$basecode";
    if (! -e $distDir) {
	print STDERR "[building]: OutDir=$distDir\n";
	make_path($distDir);
    }
    my $distFile ="$distDir/$fileName";
    my $distFile2="$distDir/$fileName.gz";
    if (-e $distFile) {
	my $size1 = -s $distFile;
	my $size2 = -s $file;

	if ($size1 == $size2) {
	    print STDERR "[Equ]: Files are same size $fileName ($size1 == $size2) skipping...\n";
	    $matCnt++;
	} else {
	    $misCnt++;
	    my $cmd="cp -f \"$file\" $distFile";
	    print STDERR "[mis]: Files are NOT same size $fileName ($size1 == $size2) skipping...\n";
	    print STDERR "[mis]: '$cmd'\n";
	    system($cmd);
	}
    } elsif (-e $distFile2) {
	my $size1 = -s $distFile2;
	my $size2 = -s $file;

	if ($size1 == $size2) {
	    print STDERR "[Equ]: Files are same size $fileName ($size1 == $size2) skipping...\n";
	    $matCnt++;
	} else {
	    $misCnt++;
	    my $cmd="cp -f \"$file\" $distFile2";
	    print STDERR "[mis]: Files are NOT same size $fileName ($size1 == $size2) skipping...\n";
	    print STDERR "[mis]: '$cmd'\n";
	    system($cmd);
	}
    } else {
	$newCnt++;
	my $cmd="cp -f \"$file\" $distFile";
	print STDERR "[new]: '$cmd'\n";
	system($cmd);

	$cmd = "gzip -f $distFile";
	system($cmd) if (!$isGZ);
    }

    # last if ($cpCnt > 1);
    # last;
}

print STDERR "[Counts]: new=$newCnt\n";
print STDERR "[Counts]: mat=$matCnt\n";
print STDERR "[Counts]: mis=$misCnt\n";
print STDERR "[Counts]: tot=$totCnt\n";

## End of file
