#! /usr/bin/perl -w

use strict;

# my $tar_file = "/Users/bretbarnes/Documents/data/workhorse/aqp_data/Epic1/auxdb/cpg-list.sorted.txt";
# my $tar_file = "/Users/bretbarnes/Documents/data/manifests/methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.snp-ids.csv";
# my $tar_file = "/Users/bretbarnes/Documents/data/workhorse/aqp_data/Epic2/auxdb/cpg-list.sorted.txt";

# Evonik v1 5k-add-on Chicken::
# my $tar_file = "/Users/bretbarnes/Documents/Projects/Envonik/add_on_5k/selected/Illumina_submission_cgids_chicken_infinium_design.txt";
#  IDX:: [19]
#  CMD:: gzip -dc Projects/Envonik/cgnDB/gal.improbe-designable.tsv.gz| /Users/bretbarnes/Documents/tools/Workhorse-Unstained/scripts/perl/bsc/subset/subset.pl > Projects/Envonik/add_on_5k/selected/extracted-designs/gal-extracted-designs.tsv
#

# Evonik v1 5k-add-on Human:
# my $tar_file = "/Users/bretbarnes/Documents/Projects/Envonik/add_on_5k/selected/Illumina_submission_cgids_human.txt";
my $tar_file = "/Users/bretbarnes/Documents/Projects/Envonik/add_on_5k/selected/Illumina_submission_cgids_human_infinium_design.txt";
#  CHange delimiter for target to COM(,) [3]; add on Infinium Design::
#  CMD:: gzip -dc /Users/bretbarnes/Documents/Projects/Envonik/add_on_5k/probes_bdf_talos_1.all_designs.23082022.v2.csv.gz |  /Users/bretbarnes/Documents/tools/Workhorse-Unstained/scripts/perl/bsc/subset/subset.pl > /Users/bretbarnes/Documents/Projects/Envonik/add_on_5k/selected/extracted-designs/hum-extracted-designs.tsv
#
#

my %tar = ();
open(IN, "<$tar_file") or die("Failed to open $tar_file for reading: $!");
while(<IN>) {
    s/^\s+//;
    s/\s+$//;

    my $line = $_;
    my @data = split(/\t/, $line);
    # my @data = split(/,/, $line);

    my $id = $data[0];

    $data[1] = 1 if ( !defined($data[1] ) );

    $tar{$id} = $data[1];;

    # print "$id\t".$tar{$id}."\n";
}
close(IN);

print STDERR "[Target]: total keys = ". scalar( keys %tar ) . "\n";

my $matCnt = 0;
my $totCnt = 0;
my $idsIdx = 0;

# Temp Fix::
$idsIdx = 2;
$idsIdx = 0;
$idsIdx = 19;
$idsIdx = 3;

while(<STDIN>) {
    s/^\s+//;
    s/\s+$//;

    my $line = $_;
    # my @data = split(/\t/, $line);
    my @data = split(/,/, $line);

    # if ($line =~ m/^Seq_ID/ || 
    if ($line =~ m/^Chromosome/ ||
	$line =~ m/^\#/ ) {
	print "$line\n";
	next;
    }

    my $id = $data[$idsIdx];

    if ( defined($tar{$id}) ) {
	print "$line\t$tar{$id}\n";

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
