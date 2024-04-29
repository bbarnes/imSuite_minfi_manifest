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

my ($id, $seq, $outDir);

my $bowExe = "/illumina/thirdparty/bowtie2/bowtie2-2.2.2/bowtie2-build";

GetOptions("o=s"  => \$outDir,

    );

make_path($outDir);


my %seqs = ();
my $lineCnt=0;

while(<STDIN>) {
    s/^\s+//;
    s/\s+$//;

    my $line = $_;

    if ($line =~ m/^>([^\s]+) .*$/) {
	$id=$1;
	# $id =~ s/\s+/_/gi;
	# $id =~ s/\|/_/gi;
	# $id =~ s/,+/_/gi;

	my $seqCnt = scalar(keys %seqs);

	$id = join("_", $id,$seqCnt);

	print STDERR "\t[Current]: ID='$id'\n";

	$seqs{$id} = "";
    } else {
	die("\nERROR: id is NULL; line($lineCnt)=$line\n") if (! defined($id));

	$seqs{$id} .= uc($line);
    }

    $lineCnt++;
}

print STDERR "Seq Count=".scalar(keys %seqs) . ", lineCount=$lineCnt\n";

# Write individual fastas::

my $outCnt=0;
foreach my $id (sort keys %seqs) {

    my $fasFH;
    my $fasFN = "$outDir/$id.fa";

    my $fasSeq = $seqs{$id};
    my $fasLen = length($fasSeq);

    if ($fasLen < 10000) {
	print STDERR "[Skipping]: id=$id; fasLen=$fasLen\n";
	next;
    }

    # open($fasFH, ">:gzip", $fasFN) or die("Failed to open $fasFN for stream gzip writing: $!");
    open($fasFH, ">$fasFN") or die("Failed to open $fasFN for stream gzip writing: $!");
    print {$fasFH} ">$id\n$seqs{$id}\n";
    close($fasFH);

    # Need to build bowtie index as well::
    #
    my $cmd="";
    if (-e $bowExe) {
	$cmd = "$bowExe $fasFN $fasFN";
	system($cmd);
    }
    $cmd = "gzip $fasFN";
    system($cmd);

    print STDERR "[Wrote]: fasFN=$fasFN\n";

    # last if ($outCnt>5);
    $outCnt++;
}

# End of file
