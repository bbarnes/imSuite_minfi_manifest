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

my @vcfFiles = ();
my ($out, $chrSplit);

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

GetOptions("v=s@" => \@vcfFiles,
	   "o=s"  => \$out,
	   "s"    => \$chrSplit,
	   );

if (! defined($out) || scalar(@vcfFiles) == 0 ) {
    print STDERR "[Usage]: $0 -v vcfFile -o out [-s] [options]\n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n";
    exit(1);
}

my $wroteFirstHeader = FALSE;
my %fhs = ();
if ($chrSplit) {
    $out =~ s/\/+$//;
    die("[ERROR]: $out is not a direcotry while using chromsome split!") if (! -d $out);
} else {
    open(OUT, ">$out")
	or die("Failed to open $out for writing!");
}
my %indexes = ();
my $totCnt  = 0;

foreach my $vcfFile (@vcfFiles) {
    open(IN, "gunzip -fc $vcfFile |")
	or die("Failed to open $vcfFile for reading: $!");

    my %svCnt = ();
    while(<IN>) {
	chomp();
	
	my $line = $_;
	if ($line =~ m/^\#/) {
	    if (! $wroteFirstHeader ) {

	    }
	    next;
	}

	my @fields = split(/\t/, $line);
	
	my $chr   = $fields[VCF_CHROM];
	my $start = $fields[VCF_POS];
	my $ref   = $fields[VCF_REF];
	my $alt   = $fields[VCF_ALT];
	my $info  = $fields[VCF_INFO];
	
	die("[ERROR]: Info field is break!") if (! defined($info) || length($info) == 0);
	
	my ($end, $svType, $impercise, $bkpt, $svLen);
	
	my @infoDat = split(/;/, $info);
	foreach my $str (@infoDat) {
	    my ($key, $val) = split(/=/, $str);
	    
	    if ($key eq "END" && $val =~ m/^[0-9-]+/) {
		$end = $val;
	    }
	    if ($key eq "SVLEN" && $val =~ m/^[0-9-]+/) {
		$svLen = $val;
	    }
	    if ($key eq "SVTYPE") {
		$svType = $val;
	    }
	    if ($key eq "BKPT") {
		$bkpt = $val;
	    }
	    if ($key eq "IMPRECISE") {
		$impercise = TRUE;
	    }
	}
	if (! defined($svType) && ! defined($bkpt)) {
	    die("[ERROR]: Info field does not contain SVTYPE or BKPT field: $info!") ;
	} elsif (! defined($svType)) {
	    die("[ERROR]: BKPT defined with END: $info!") if (defined($end));
	    $end = $start + 1;
	    $svType = "BKPT:$bkpt";
	}
	die("[ERROR]: Info field does not contain SVTYPE field: $info!") if (! defined($svType));
	die("[ERROR]: BKPT and both IMPERCISE not defined: $info!") 
	    if ((! defined($bkpt) && defined($impercise)) ||
		(! defined($impercise) && defined($bkpt)));
	
	die("[ERROR]: Both END and SVLEN are not defined: $info") if (! defined($end) && ! defined($svLen));
	die("[ERROR]: START, END and SVLEN do not agree: $start, $end, $svLen: $line!") 
	    if (defined($svLen) && defined($end) && $start + $svLen - 1 != $end);
	if (defined($svLen)) {
	    die("[ERROR]: Failed logic! $line") if (defined($end));
	    if ($svType eq "INS") {
		$end = $start + 1;
	    } elsif ($svType eq "BND") {
		$end = $start + 1;
	    } else {
		$end = $start + abs($svLen) - 1;
	    }
	}
	next if ($line =~ m/GROUPERuTRAN/);
	die("[ERROR]: Info field does not contain END field: $info!") if (! defined($end));
	die("[ERROR]: Confusing coordinates for INS: $line!") if ($svType eq "INS" && $start != $end - 1);
	die("[ERROR]: START > END ($start > $end): $line") if ($start > $end);
	
	next if (! defined($knownTypes{$svType}));
	
	print STDERR "[Warning]: Unknown type: $svType, on line: $line\n"  if (! defined($knownTypes{$svType}));
	
	$end++ if ($start == $end);
	
	if (! defined($svCnt{$svType})) {
	    $svCnt{$svType} = 1;
	} else {
	    $svCnt{$svType}++;
	}
	
	## Force insertion length to be added
	if ($svType eq "INS") {
	    my $insLen = length($alt) - length($ref);
	    die("[ERROR]: Failed matching insertion length: $insLen != $svLen: line: $line!") 
		if (defined($svLen) && $svLen != $insLen);
	    $svType .= ":$svLen";
	}

	if ($chrSplit) {
	    if (! defined($fhs{$chr})) {
		my $outFile = "$out/$chr.bed";
		open($fhs{$chr}, ">$outFile")
		    or die("Failed to open $outFile for writing: $!");
	    }
	    print {$fhs{$chr}} "$chr\t$start\t$end\t$idx\t$svType\n";
	} else {
	    print OUT "$chr\t$start\t$end\t$idx\t$svType\n";
	}
	$idx++;
    }
    close(IN);
    $wroteFirstHeader = TRUE;
}

if ($chrSplit) {
    foreach my $key (keys %fhs) {
	close($fhs{$key});
    }
} else {
    close(OUT);
}

print STDERR "\# Processed $idx variants\n";

foreach my $key (sort keys %svCnt) {
    print STDERR "\# $key\t$svCnt{$key}\n";
}

## End of file
