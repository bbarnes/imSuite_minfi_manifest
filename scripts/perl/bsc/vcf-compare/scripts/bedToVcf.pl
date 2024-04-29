#! /usr/bin/perl -w

use strict;
use Getopt::Long;

use constant FALSE => 0;
use constant TRUE  => 1;

use constant BED_CHROM  => 0;
use constant BED_START  => 1;
use constant BED_END    => 2;
use constant BED_STRAND => 3;
use constant BED_INFO   => 4;

use constant VCF_CHROM  => 0;
use constant VCF_POS    => 1;
use constant VCF_ID     => 2;
use constant VCF_REF    => 3;
use constant VCF_ALT    => 4;
use constant VCF_QUAL   => 5;
use constant VCF_FILTER => 6;
use constant VCF_INFO   => 7;

my @bedFiles = ();
my ($headerFile);

$headerFile = "/bioinfoSD/users/bbarnes/Io/CNV/headers/1000Genomes.header.txt";

my %knownTypes = ("DEL"        => 1,
		  "INS"        => 1,
		  "INV"        => 1,
		  "DUP"        => 1,
		  "DUP:TANDEM" => 1,
		  "ADJ"        => 1,
		  "CNV"        => 1,
		  "BND"        => 1,
		  "BKPT"       => 1,
		  "BKPT:LEFT"  => 1,
		  "BKPT:RIGHT" => 1,
		  );

GetOptions("b=s@" => \@bedFiles,
	   "h=s"  => \$headerFile,
	   );

if (! defined($headerFile) || scalar(@bedFiles) == 0) {
    print STDERR "[Usage]: $0 -b bedFile -h headerFile [options]\n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n";
    exit(1);
}

my @header = ();

print STDERR "[Loading]: headerFile=$headerFile\n";
open(IN, "<$headerFile")
    or die("Failed to open $headerFile for reading: $!");
while(<IN>) {
    s/^\s+//;
    s/\s+$//;
    next if /^(\s)*$/;
    #next if /^#/;

    my $line = $_;

    print "$line\n";
}
close(IN);
print STDERR "[Header]: Wrote header!\n";

my $totCnt  = 0;
foreach my $bedFile (@bedFiles) {
    open(IN, "gunzip -fc $bedFile |")
	or die("Failed to open $bedFile for reading: $!");

    my %svCnt = ();
    while(<IN>) {
	chomp();
	
	my $line = $_;
	my @data = split(/\t/, $line);

	my @vcf = (); 
	my ($vcfEnd, $ref, $alt, $qual, $filt, $info, $format, $sample);

	$vcf[VCF_CHROM]  = $data[BED_CHROM];
	$vcf[VCF_POS]    = $data[BED_START];
	$vcfEnd          = $data[BED_END];
	#$vcf[VCF_STRAND] = $data[BED_STRAND];
	$vcf[VCF_ID]     = $data[BED_INFO];

	$ref  = "N";
	$alt  = "<DEL>";
	$qual = "1000";
	$filt = "PASS";
	$info = "END=$vcfEnd;SVTYPE=DEL";
	$format = "GT";
	$sample = "0/1";

	print "$vcf[VCF_CHROM]\t$vcf[VCF_POS]\t$vcf[VCF_ID]\t$ref\t$alt\t$qual\t$filt\t$info\t$format\t$sample\n";
	$totCnt++;
    }
}
close(IN);

print STDERR "[Done]: totalCount=$totCnt\n";

## End of file
