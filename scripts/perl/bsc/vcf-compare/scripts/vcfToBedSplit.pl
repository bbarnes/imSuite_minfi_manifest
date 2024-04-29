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
use constant VCF_SAMPLE_FMT => 8;
use constant VCF_SAMPLE_VAL => 9;

my ($out, $chrSplit);

my $unkTypeCnt = 0;
my $trnSkipCnt = 0;

my $delSizeLimit = 0;

my $cnvOnly;

my %knownTypes = ("DEL"             => 1,
		  "INS"             => 1,
		  "INV"             => 1,
		  "DUP"             => 1,
		  "DUP:TANDEM"      => 1,
		  "TRN"             => 1,
		  "CNV"             => 1,
		  "BND:OPENEND"     => 1,
		  "BND:OPENEND:DEL" => 1,
		  "BND:OPENEND:INS" => 1,
		  "BND:OPENEND:INV" => 1,
		  "BND:OPENEND:DUP:TANDEM" => 1,
		  "BND:OPENEND:TRN" => 1,
		  "COMPLEX"         => 1,
		  "BND"             => 1,
		  );

my %svCounts = ();
foreach my $key (keys %knownTypes) {
    $svCounts{$key} = 0;
}

GetOptions("o=s"  => \$out,
	   "s"    => \$chrSplit,
	   "d=i"  => \$delSizeLimit,
	   "cnvOnly" => \$cnvOnly,
	   );

if (! defined($out) || ! defined($delSizeLimit)) {
    print STDERR "[Usage]: $0 -o out -d deletionSizeLimit [-s] [-cnvOnly] [options] < *.vcf \n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n";
    exit(1);
}

print STDERR "[Using] cnvOnly=TRUE!\n" if (defined($cnvOnly));

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

my $bkptCnt = 0;
while(<STDIN>) {
    s/^\s+//;
    s/\s+$//;
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
    @fields = @fields[0...7];
	
    my $chr   = $fields[VCF_CHROM];
    my $start = $fields[VCF_POS];
    my $vcfId = $fields[VCF_ID];
    my $ref   = $fields[VCF_REF];
    my $alt   = $fields[VCF_ALT];
    my $info  = $fields[VCF_INFO];
    my $forVal = $fields[VCF_SAMPLE_VAL];

    #next if ($forVal !~ m/:/);
	
    die("[ERROR]: Info field is break! Line: $line") if (! defined($info) || length($info) == 0);
	
    my ($end, $svType, $svLen);
	
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
    }

    next if (! defined($svType) || length($svType) == 0);
    next if (defined($cnvOnly) &&
	     ($svType ne "DEL" || $svType =~ m/^DUP/ || $svType =~ m/^CNV/));

    ## Yet again another special fix for skipping SNVs from simulation that have
    ##  no SVTYPE
    #next if (length($ref) == length($alt) && $vcfId =~ m/_rs/);
    next if ($vcfId =~ m/_rs[0-9]+$/);

    ## This is for skipping SNPs 
    next if (! defined($svType) && length($ref) == length($alt));
    next if ($alt =~ m/,/ || $ref =~ m/,/);

    ## Fix for missing SVTYPE for Starling calls
    if (! defined($svType) && $info =~ m/CIGAR=[^;]*D[^;]*;/) {
	$svType = "DEL";
	$svLen  = length($alt) - length($ref);
	$end    = $start + abs($svLen);
    }
    if (! defined($svType) && $info =~ m/CIGAR=[^;]*I[^;]*;/) {
	$svType = "INS";
	$svLen  = length($alt) - length($ref);
	$end    = $start + abs($svLen);
    }
    ##
    ## FIXME: Temp fix for canvas data
    ##
    next if (! defined($svLen));

    ## Require an SVTPE to be defined!
    die("[ERROR]: No SV TYPE defined for variant: $line!") if (! defined($svType));

    ## Skip deletions that are small and overlap a large region.
    ##  This makes picking rc values easier...
    next if ($svType eq "DEL" && abs($svLen) < $delSizeLimit);
    ## Actually we'll do it for everything...
    next if ($svType ne "BND" && abs($svLen) < $delSizeLimit);

    ## First fix the BND types
    if ($svType =~ m/BND/) {
	if ($alt =~ m/\[/ || $alt =~ m/\]/) {
	    $svType = "TRN";
	    $svLen  = 1;
	} else {
	    $svType = $alt;
	    $svType =~ s/<//gi;
	    $svType =~ s/>//gi;
	}
    }

    ## Temp fix for Isis
    if ($svType =~ m/INV:TANDEM/) {
	$svType =~ s/:TANDEM//;
    }

    ## Added TANEM to duplications
    $svType .= ":TANDEM" if ($svType eq "DUP");

    ## Ensure translocations are annotated correctly
    if ($vcfId =~ m/GROUPERuTRAN/) {
	die("Grouper type is translocation, but ALT and SVTYPE are not on line: $line\n")
	    if ($svType ne "TRN");
    }

    ## Skip SNVs
    next if ($svType eq "SNV");

    ## Complain if SV type is unknown
    if (! defined($knownTypes{$svType})) {
	print STDERR "[Warning]: Unknown SV Type: $svType, on line: $line! Skipping...\n";
	if (! defined($svCounts{"UNK"})) {
	    $svCounts{"UNK"} = 1;
	} else {
	    $svCounts{"UNK"}++;
	}
	next;
    }
    $svCounts{$svType}++;

    ## CNV size calc
    if ($svType eq "CNV") {
	$svLen = abs($end - $start);
	#$line =~ s/SVTYPE=CNV/SVTYPE=CNV;SVLEN=$svLen/;
    }

    ## Special handleing for new inversion code
    ##  NOT SURE WHAT THIS DOES ANYMORE...
    #if (defined($eventType) && $eventType =~ m/^INV([0-9]+)$/) {
	#next if (defined($svIds{$eventType}));
	#$svIds{$eventType} = 1;
	#$svType = "INV";
	#$line =~ s/SVTYPE=BND/SVTYPE=INV/;
    #}
    
    die("[ERROR]: Both END and SVLEN are not defined: $info\nLine: $line!") if (! defined($end) && ! defined($svLen));
    my $calcLen = $start + abs($svLen);
    #my $diffInLen = $calcLen - $end;
    ## FIXME: Tandem Dups in Eagle do not have the correct start and ened coordinates. There is a fix here to make that simply pass through
    die("[ERROR]: START, END and SVLEN do not agree: $start + $svLen ($calcLen) != $end:\nLine: $line!") 
	if (defined($svLen) && defined($end) && $calcLen != $end 
	    && $svType ne "INS" && $svType =~ m/BND/ && $svType ne "TRN" && $svType ne "DUP:TANDEM");
    if (defined($svLen)) {
	if ($svType eq "INS" || $svType =~ m/BND/ || $svType eq "TRN") {
	    $end = $start + 1;
	} else {
	    $end = $start + abs($svLen);
	}
    }
    die("[ERROR]: Info field does not contain END field: $info!") if (! defined($end));
    #die("[ERROR]: Confusing coordinates for INS: $line!") if (($svType eq "INS" || $svType =~ m/BND/ || $svType eq "TRN") && $start != $end - 1);
    die("[ERROR]: START > END ($start > $end): $line") if ($start > $end);    
    
    ## This should never happen now...
    $end++ if ($start == $end);
	
    ## Force insertion length to be added
    my $svTypeLen = $svType;
    if ($svType eq "INS") {
	my $insLen = length($alt) - length($ref);
	if (defined($svLen) && $svLen != $insLen && $svType ne "BKPT") {
	    #print STDERR "[Warning]: Failed matching insertion length: $insLen != $svLen: line: $line!\n";
	    next;
	}
	#die("[ERROR]: Failed matching insertion length: $insLen != $svLen: line: $line!") 
	#    if (defined($svLen) && $svLen != $insLen && $svType ne "BKPT");
	$svTypeLen .= ":$svLen";
    }
    
    if ($chrSplit) {
	## Leaving this code for now, but it should not be used
	die("This option (split by chromosome) is no longer supported!");
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
	$fields[VCF_INFO] =~ s/;$//;
	$fields[VCF_INFO] .= ";SVTYPE_C=$svType;SVLEN_C=$svLen";
	$line = join("\t", @fields);

	print {$fhsVcf{$chr}} "$line\n";
	print {$fhsBed{$chr}} "$chr\t$start\t$end\t$cntVcf{$chr}\t$svTypeLen\n";
    } else {
	$fields[VCF_INFO] =~ s/;$//;
	$fields[VCF_INFO] .= ";SVTYPE_C=$svType;SVLEN_C=$svLen";
	$line = join("\t", @fields);

	print OUTVCF "$line\n";
	print OUTBED "$chr\t$start\t$end\t$totIdx\t$svTypeLen\n";
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

print STDERR "\# Processed $totIdx variants\n";
foreach my $key (sort keys %svCounts) {
    print STDERR "\# $key\t$svCounts{$key}\n";
}
print STDERR "\n";

## End of file
