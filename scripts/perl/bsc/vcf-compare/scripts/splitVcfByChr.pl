#! /usr/bin/perl -w

use strict;
use Getopt::Long;

use constant FALSE => 0;
use constant TRUE  => 1;

use constant VCF_CHROM  => 0;
use constant VCF_POS    => 1;
use constant VCF_ID     => 2;
use constant VCF_REF    => 3;
use constant VCF_ALT    => 4;
use constant VCF_QUAL   => 5;
use constant VCF_FILTER => 6;
use constant VCF_INFO   => 7;

my ($vcfFile, $outDir);

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

GetOptions("v=s"  => \$vcfFile,
	   "o=s"  => \$outDir,
	   );

if (! defined($vcfFile)) {
    print STDERR "[Usage]: $0 -v vcfFile -o out [options]\n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n";
    exit(1);
}

open(IN, "gunzip -fc $vcfFile |")
    or die("Failed to open $vcfFile for reading: $!");

my %fhs = ();
$outDir =~ s/\/+$//;
die("[ERROR]: $outDir is not a direcotry while using chromsome split!") if (! -d $outDir);

my $header = "";
my $idx;
while(<IN>) {
    chomp();

    my $line = $_;
    if ($line =~ m/^\#/) {
	$header .= $line . "\n";
	next;
    }
    my @fields = split(/\t/, $line);

    my $chr   = $fields[VCF_CHROM];
    if (! defined($fhs{$chr})) {
	my $outFile = "$outDir/$chr.vcf";
	open($fhs{$chr}, ">$outFile")
	    or die("Failed to open $outFile for writing: $!");
	print {$fhs{$chr}} $header;
    }
    print {$fhs{$chr}} "$line\n";
    $idx++;
}
close(IN);

foreach my $key (keys %fhs) {
    close($fhs{$key});
}

print STDERR "\# Processed $idx variants\n";

## End of file
