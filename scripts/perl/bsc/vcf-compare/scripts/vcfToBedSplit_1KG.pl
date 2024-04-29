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

my ($out, $chrSplit);

my $unkTypeCnt = 0;
my $trnSkipCnt = 0;

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

GetOptions(
	   "o=s"  => \$out,
	   "s"    => \$chrSplit,
	   );

if (! defined($out)) {
    print STDERR "[Usage]: $0 -o out [-s] [options] < *.vcf \n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n";
    exit(1);
}

my $loadedHeader = FALSE;
my %cntVcf = ();
my %fhsBed = ();
my %fhsVcf = ();
if ($chrSplit) {
    $out =~ s/\/+$//;
    die("[ERROR]: $out is not a direcotry while using chromsome split!") if (! -d $out);
} else {
    open(OUTBED, ">$out.bed")
	or die("Failed to open $out.bed for writing!");
    open(OUTVCF, ">$out.vcf")
	or die("Failed to open $out.vcf for writing!");
}
my %indexes = ();
my $totIdx  = 0;
my $header = "";

my %svCnt = ();
my %svIds = ();
while(<STDIN>) {
    chomp();
	
    my $line = $_;
    if ($line =~ m/^\#/) {
	if (! $loadedHeader ) {
	    $header .= $line . "\n";
	    print OUTVCF $line . "\n" if (! $chrSplit)
	}
	next;
    }
    $loadedHeader = TRUE;

    my @fields = split(/\t/, $line);
	
    my $chr   = $fields[VCF_CHROM];
    my $start = $fields[VCF_POS];
    my $ref   = $fields[VCF_REF];
    my $alt   = $fields[VCF_ALT];
    my $info  = $fields[VCF_INFO];
	
    die("[ERROR]: Info field is break! Line: $line") if (! defined($info) || length($info) == 0);
	
    my ($end, $svType, $impercise, $bkpt, $svLen, $eventType);
	
    my @infoDat = split(/;/, $info);
    foreach my $str (@infoDat) {
	my ($key, $val) = split(/=/, $str);
	
	if ($key eq "END" && $val =~ m/^[0-9-]+/) {
	    $end = $val;
	}
	if ($key eq "SVLEN" && $val =~ m/^[-0-9-]+/) {
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
	if ($key eq "EVENT") {
	    $eventType = $val;
	}
    }

    ## Skip SNVs
    next if ($svType eq "SNV");

    if (! defined($svType) && ! defined($bkpt)) {
	die("[ERROR]: Info field does not contain SVTYPE or BKPT field: $info!") ;
    } elsif (! defined($svType)) {
	die("[ERROR]: BKPT defined with END: $info!") if (defined($end));
	$end = $start + 1;
	$svType = "BKPT:$bkpt";
    }
    die("[ERROR]: Info field does not contain SVTYPE field: $info!") if (! defined($svType));
    #die("[ERROR]: BKPT and both IMPERCISE not defined: $info!\nLine: $line\n") 
    #if ((! defined($bkpt) && defined($impercise)) ||
    #(! defined($impercise) && defined($bkpt)));

    if ((! defined($bkpt) && defined($impercise)) ||
	(! defined($impercise) && defined($bkpt))) {
	print STDERR ("[ERROR]: BKPT and both IMPERCISE not defined: $info!\nLine: $line\n");
	    next;
    }
    
    ## Skip translocations
    if ($line =~ m/GROUPERuTRAN/ ||
	(defined($eventType) && $eventType =~ m/^RR/)) {
	$trnSkipCnt++;
	next;
    }

    ## Special handleing for new inversion code
    if (defined($eventType) && $eventType =~ m/^INV([0-9]+)$/) {
	next if (defined($svIds{$eventType}));
	$svIds{$eventType} = 1;
	$svType = "INV";
	$line =~ s/SVTYPE=BND/SVTYPE=INV/;
    }
    
    ## Swap start and end for special 1000G case
    #my $tmpPos = $start;
    #$start = $end;
    #$end = $tmpPos;
    #$svLen = abs($svLen);

    die("[ERROR]: Both END and SVLEN are not defined: $info\nLine: $line!") if (! defined($end) && ! defined($svLen));
    my $calcLen = $start + abs($svLen);
    die("[ERROR]: START, END and SVLEN do not agree: $start + $svLen ($calcLen) != $end:\nLine: $line!") 
	if (defined($svLen) && defined($end) && $calcLen != $end && $svType ne "INS" && $svType ne "BND");
    if (defined($svLen)) {
	#die("[ERROR]: Failed logic! $line") if (defined($end));
	if ($svType eq "INS") {
	    $end = $start + 1;
	} elsif ($svType eq "BND") {
	    $end = $start + 1;
	} else {
	    $end = $start + abs($svLen);
	}
    }
    die("[ERROR]: Info field does not contain END field: $info!") if (! defined($end));
    die("[ERROR]: Confusing coordinates for INS: $line!") if ($svType eq "INS" && $start != $end - 1);
    die("[ERROR]: START > END ($start > $end): $line") if ($start > $end);
    
    if (! defined($knownTypes{$svType})) {
	$unkTypeCnt++;
	next;
    }
    
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
	if (! defined($fhsBed{$chr})) {
	    my $outFile = "$out/$chr.bed";
	    my $vcfFile = "$out/$chr.vcf";
	    print STDERR "\# Opening output files for $chr:\n";
	    print STDERR "\#\t$outFile\n";
	    print STDERR "\#\t$vcfFile\n";
	    open($fhsBed{$chr}, ">$outFile")
		or die("Failed to open $outFile for writing: $!");
	    open($fhsVcf{$chr}, ">$vcfFile")
		or die("Failed to open $vcfFile for writing: $!");
	    print {$fhsVcf{$chr}} $header;
	    $cntVcf{$chr} = 0;
	}
	## Set counting index
	print {$fhsVcf{$chr}} "$line\n";
	print {$fhsBed{$chr}} "$chr\t$start\t$end\t$cntVcf{$chr}\t$svType\n";
    } else {
	print OUTVCF "$line\n";
	print OUTBED "$chr\t$start\t$end\t$totIdx\t$svType\n";
    }
    $totIdx++;
    $cntVcf{$chr}++;
}

if ($chrSplit) {
    foreach my $chr (keys %fhsBed) {
	print STDERR "\# Loaded: " . $cntVcf{$chr} . " variants for chr: $chr\n";
	close($fhsBed{$chr});
	close($fhsVcf{$chr});
    }
} else {
    close(OUTBED);
    close(OUTVCF);
}

print STDERR "\# Processed $totIdx variants, skipped $unkTypeCnt\n";
print STDERR "\#  Skipped $trnSkipCnt translocations for now...\n";

foreach my $key (sort keys %svCnt) {
    print STDERR "\# $key\t$svCnt{$key}\n";
}

## End of file
