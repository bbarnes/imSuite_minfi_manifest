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
my $type="idat";

my %dat=();

GetOptions("s=s" => \$srcDir,
	   "o=s" => \$outDir,
	   "t=s" => \$type,
    );

if (! defined($srcDir) || ! defined($outDir)) {
    print STDERR "[Usage]: $0 -t type -s srcDir -o outDir < find-list [options]\n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n";
    exit(1);
}

die("[ERROR]: Source Directory must exist='$srcDir'\n\n") if (! -e $srcDir);
$srcDir =~ s/\/+$//;
$outDir =~ s/\/+$//;
if (! -e $outDir) { make_path($outDir); }

my $cpCnt = 0;
while(<STDIN>) {
    s/^\s+//;
    s/\s+$//;

    my $file=$_;
    $file =~ s/^\.\///;

    $file = $srcDir."/".$file;
    die("[ERROR]: file=$file does not exist!\n") if (! -e $file);

    my ($basecode, $poscode, $fileName, $color);

    if ($type eq "idat") {
	if ($file =~ m/^.*\/(([0-9]+)_([RC0-9]+)_(Grn|Red)\.idat.gz)/) {
	    $fileName=$1;
	    $basecode=$2;
	    $poscode=$3;
	    $color=$4;
	} else {
	    print STDERR "[Skipping]: Unrecognized IDAT file=$file\n";
	}
    } else {
	if ($file =~ m/^.*\/(([0-9]+)_([RC0-9]+)_([0-9]+)\.dmap.gz)/) {
	    $fileName=$1;
	    $basecode=$2;
	    $poscode=$3;
	    $color=$4;
	} else {
	    print STDERR "[Skipping]: Unrecognized DMAP file=$file\n";
	}
    }
    if (defined($dat{$fileName})) {
	print STDERR "[Duplicate]: $fileName\n";
	next;
    }

    print STDERR "fileName=$file, basecode=$basecode, poscode=$poscode\n";

    my $distDir = $outDir."/$basecode";
    if (! -e $distDir) {
	print STDERR "[building]: OutDir=$distDir\n";
	make_path($distDir);
    }
    my $cmd="cp -f $file $distDir/";
    print STDERR "[cmd]: '$cmd'\n";
    system($cmd);


    $cpCnt++;
    #last if ($cpCnt > 1);
}


## End of file
