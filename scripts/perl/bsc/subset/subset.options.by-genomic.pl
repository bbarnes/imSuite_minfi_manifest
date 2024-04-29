#! /usr/bin/perl -w

use strict;
use Getopt::Long;

# my $tar_file = "/Users/bretbarnes/Documents/data/workhorse/aqp_data/Epic1/auxdb/cpg-list.sorted.txt";
# my $tar_file = "/Users/bretbarnes/Documents/data/manifests/methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.snp-ids.csv";
# my $tar_file = "/Users/bretbarnes/Documents/data/workhorse/aqp_data/Epic2/auxdb/cpg-list.sorted.txt";

my $tar_file;

GetOptions(
    "t=s" => \$tar_file,
    );

if (! defined($tar_file) ) {
    print STDERR "[Usage]: $0 -t target_file.txt < STDIN > STDOUT\n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n\n";
    exit(1);
}

my %tar = ();
open(IN, "<$tar_file") or die("Failed to open $tar_file for reading: $!");
while(<IN>) {
    s/^\s+//;
    s/\s+$//;

    my $line = $_;
    my @data = split(/\t/, $line);

    my $id = $data[0];
    $id =~ s/^chr//;

    $tar{$id} = 1;
}
close(IN);

print STDERR "[Target]: total keys = ". scalar( keys %tar ) . "\n";

my $matCnt = 0;
my $totCnt = 0;
my $idsIdx = 0;

$idsIdx = 2;
$idsIdx = 0;

while(<STDIN>) {
    s/^\s+//;
    s/\s+$//;

    my $line = $_;
    my @data = split(/\t/, $line);

    if ($line =~ m/^Seq_ID/ ||
	$line =~ m/^\#/ ) {
	print "$line\n";
	next;
    }

    my $id = $data[$idsIdx];
    my $pos = "$data[3],$data[4]";
    
    if ( defined($tar{$pos}) ) {
	print "$line\n";

	$matCnt++;
    }
    if ($totCnt % 1000000 == 0) {
	print STDERR "[Status]: $matCnt, $totCnt\n";
    }
    
    $totCnt++;
}
close(IN);

print STDERR "[Done]: $matCnt/$totCnt\n";

# End of file
