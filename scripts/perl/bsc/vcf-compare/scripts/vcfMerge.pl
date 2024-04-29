#!/usr/bin/env perl

# Copyright (c) Illumina 2010
# Author: Bret Barnes & Han-Yu Chuang
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.

use warnings;
use strict;

use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;

use constant FALSE => 0;
use constant TRUE  => 1;

## VCF File Formats
use constant VCF_CHROM   => 0;
use constant VCF_POS     => 1;
use constant VCF_ID      => 2;
use constant VCF_REF     => 3;
use constant VCF_ALT     => 4;
use constant VCF_QUAL    => 5;
use constant VCF_FILTER  => 6;
use constant VCF_INFO    => 7;
use constant VCF_FORMAT  => 8;
use constant VCF_SAMPLE  => 9;

## Intersection File Formats
use constant INT_CHR1   => 0;
use constant INT_START1 => 1;
use constant INT_END1   => 2;
use constant INT_IDX1   => 3;
use constant INT_TYPE1  => 4;
use constant INT_CHR2   => 5;
use constant INT_START2 => 6;
use constant INT_END2   => 7;
use constant INT_IDX2   => 8;
use constant INT_TYPE2  => 9;

## BED File formats
use constant BED_CHR   => 0;
use constant BED_START => 1;
use constant BED_END   => 2;
use constant BED_IDX   => 3;
use constant BED_INFO  => 4;
use constant BED_CALLS => 5;

## Failure codes
use constant PASS_CALL    => 1;
use constant MISRP_CALL   => 2;
use constant BKPT_CALL    => 3;
use constant MISTYPE_CALL => 4;
use constant WORST_CALL   => 5;

## Binning params
use constant BIN_START => 0;
use constant BIN_END   => 12;

## Global variables
my %options = ();
my %params  = ();

$options{'rcOverlap'}  = 0.90;
$options{'minDist'}    = 10;
$options{'seqBuf'}     = 6;

$options{'verbose'} = TRUE;
$options{'vcfLenThresh'} = 50;

my $pAll = TRUE;

my %allData = ();

## Variable types 
my @ranges = ("1,9",
              "10,49",
              "50,249",
              "250,999",
              "1000,2499",
              "2500,9999",
	      "10000,19999",
	      "20000,49999",
	      "50000,100000"
              );

#@ranges = ("1,299",
#	   "300,20000"
#	   );

my $lastBinVal = $ranges[scalar(@ranges) - 1];
$lastBinVal =~ s/^.*:([0-9]+)$/$1/;

my @svTypes = ("DEL",
	       "INS",
	       "INV",
	       "DUP:TANDEM",
	       "CNV",
	       "TRN",
	       "BND:OPENEND:DEL",
	       "BND:OPENEND:DUP:TANDEM",
	       "BND:OPENEND:INS",
	       "BND:OPENEND:INV",
	       "BND:OPENEND:TRN",
	       "COMPLEX"
	       );
#	       "TRN");

my @svCalls = (1, 2, 3, 4, 5);

## Build hash lookup
my %hist = ();
my %map  = ();
for (my $gCnt = 0; $gCnt < @ranges; $gCnt++) {
    my @a = split(/,/, $ranges[$gCnt]);
    my $start = $a[0];
    my $end   = $a[1];

    for (my $ii = $start; $ii <= $end; $ii++) {
        $map{$ii} = $gCnt;
    }
    for (my $ii = BIN_START; $ii <= BIN_END; $ii++) {
        $hist{$gCnt}[$ii] = 0;
    }
}
my $lastBin = @ranges;
for (my $ii = BIN_START; $ii <= BIN_END; $ii++) {
    $hist{$lastBin}[$ii] = 0;
    my $openBin = $lastBin + 1;
    $hist{$openBin}[$ii] = 0;
}
push @ranges, "> $lastBinVal";

#------------------------------------------------------------------------------
#
# Program Subroutines
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

sub processArgs()
{
    my $usage =
	"Usage: $0 <options>\n"
	. "\t--can=PATH(s)              - Path(s) to candidate BED calls.\n"
	. "\t--ref=PATH(s)              - Path(s) to reference BED calls.\n"
	. "\t--int=PATH(s)              - Path(s) to intersection files. It is assumed that these match the\n"
	. "\t                              order of candidates to reference BED files on the command line\n"

	. "\t--outDir=PATH              - Output directory\n"

        . "Optional :-\n"
	. "\t--vcfLenThresh             - Minimum size theshold to consider an SV (default = $options{'vcfLenThresh'})\n"
	. "\t--seqComp                  - Do sequence comparison for insertions and deletions (default = FALSE)\n"
	. "\t--rc=FLOAT                 - Reciprocal overlap percentatge threshold (default = $options{'rcOverlap'})\n"
	. "\t--md=INT                   - Minimum Distance mapping threshold (defualt = $options{'minDist'})\n"
	. "\t--sb=INT                   - Up and downstream length to use for sequence comparison. Larger values require\n"
	. "\t                              more matching, but may cover longer homopolymer cases. A value of -1 will\n"
	. "\t                              use actual up and downstream sequence lengths (default = $options{'seqBuf'})\n"

	. "\t--debug                    - Run debug mode (default = false).\n"

	. "\t-v|--verbose               - Run with verbose on.\n"
        . "\t-h|--help                  - Prints this message.\n";

    my $optionsResult
	= GetOptions(\%options,
		     "can-files|can|c=s@",
		     "ref-files|ref|r=s@",
		     "int-files|int|i=s@",

		     "vcfLenThresh=i",

		     "out-dir|outDir|out|o=s",

		     "ref-filt-thresh|refThresh|rt=f",
		     "can-filt-thresh|canThresh|ct=f",

		     "rcOverlap|rc=f",
		     "minDist|md=i",
		     "seqBuf|sb=i",

		     "noParents|np",
		     "seqComp",

		     "debug",
		     "help|h",
		     "verbose|v",
		     );
    
    if ($options{'help'}) {
        print STDERR $usage, "\n";
        exit(0);
    }

    ## Check inputs
    die($usage) if (! defined($options{'can-files'}) || @{$options{'can-files'}} == 0 ||
		    ! defined($options{'ref-files'}) || @{$options{'ref-files'}} == 0 ||
		    ! defined($options{'int-files'}) || @{$options{'int-files'}} == 0 ||
		    ! defined($options{'out-dir'}));

    ## Just quickly check all files exists
    foreach my $file (@{$options{'can-files'}}) {
	if (! -e $file) {
	    die("Candidate file: $file does not exist!");
	}
	mssg("\# Will process candidate file: $file\n");
    }
    mssg("\# Will process a total of " . @{$options{'can-files'}} . " candidate files\n");
    mssg("\n");

    foreach my $file (@{$options{'ref-files'}}) {
	if (! -e $file) {
	    die("Reference file: $file does not exist!");
	}
	mssg("\# Will process candidate file: $file\n");
    }
    mssg("\# Will process a total of " . @{$options{'ref-files'}} . " reference files\n");

    # Build work director
    buildDir($options{'out-dir'});

    mssg("\n");
}

#------------------------------------------------------------------------------

sub debug($)
{
    my ($str) = @_;

    if ($options{'debug'}) {
        print STDERR $str;
    }
}

#------------------------------------------------------------------------------

sub mssg($)
{
    my ($str) = @_;

    if ($options{'verbose'}) {
	print STDERR $str;
    }
}

#------------------------------------------------------------------------------

sub buildDir($)
{
    my ($dir) = @_;

    if (-e $dir) {
        my $clean_cmd = "rm -rf $dir/*";
        mssg("\# Cleaning '$dir'\n");
        system($clean_cmd);
    } else {
        my @dir_array = split(/\//, $dir);
        my $tmp_dir = "";
        foreach my $val (@dir_array) {
            next if ($val =~ /^(\s)*$/);

            if (! $tmp_dir) {
                if ($dir =~ m/^\/.*$/) {
                    $tmp_dir = "/" . $val;
                } else {
                    $tmp_dir = $val;
                }
            } else {
                $tmp_dir .= "/" . $val;
            }

            if (! -e $tmp_dir) {
                my $mk_cmd = "mkdir $tmp_dir";

                mssg("\# Building '$tmp_dir'\n");
                system($mk_cmd);
            }
        }
    }
}

#------------------------------------------------------------------------------

sub writeArray($$$)
{
    my ($dat, $fh, $sep) = @_;

    foreach my $val (@{$dat}) {
	print $fh $val . $sep;
    }
}

#------------------------------------------------------------------------------

sub validSequences($$)
{
    my ($can, $ref) = @_;

    my @ret = (FALSE, "-", "-");

    return(@ret) if ($$can[BED_INFO]{'Variant_length'} != $$ref[BED_INFO]{'Variant_length'});

    my $canUp  = $$can[BED_INFO]{'Upstream'};
    my $canVar = $$can[BED_INFO]{'Variant_seq'};
    my $canDn  = $$can[BED_INFO]{'Downstream'};

    my $refUp  = $$ref[BED_INFO]{'Upstream'};
    my $refVar = $$ref[BED_INFO]{'Variant_seq'};
    my $refDn  = $$ref[BED_INFO]{'Downstream'};

    if (! defined($canVar)) { $canVar = ""; }
    if (! defined($refVar)) { $refVar = ""; }

    my $minUp = $options{'seq-buf'};
    my $minDn = $options{'seq-buf'};

    my $dif = 0;
    if (! defined($$ref[BED_START])) {
	die("No start ref defined: " . bedToStr(\@{$$ref}));
    }
    if (! defined($$can[BED_START])) {
	die("No start can defined: " . bedToStr(\@{$$can}));
    }
    $dif = $$ref[BED_START] - $$can[BED_START];

    ## Find minumum up and down stream lengths for substring
    if (length($canUp) < $minUp) {
	$minUp = length($canUp);
    }
    if (length($canDn) < $minDn) {
	$minDn = length($canDn);
    }
    if (length($refUp) < $minUp) {
	$minUp = length($refUp);
    }
    if (length($refDn) < $minDn) {
	$minDn = length($refDn);
    }

    ## Cases:
    ##  1. Closed insertion or deletion
    ##  2. Closed inversion or duplication
    ##  3. Open insertion or deletion
    ##  4. Open inversion
    ##  5. Open duplication
    ##  6. Open open

    ## Case 1: Only one implmeneted so far.. [FIXME]
    ## For complete insertion or deletion construct sequecnes
    ##  and do an exact match compariosn. 

    if (! defined($dif)) {
	die("No start both defined: " .bedToStr(\@{$$ref}). "\n" . bedToStr(\@{$$can}));
    }

    if (FALSE && ! defined(substr($canUp, length($canUp) - $minUp + $dif, $minUp - $dif)) ||
	! defined(substr($canDn, 0, $minDn + $dif))) {
	print STDERR "\# DIF: $dif, $minUp, $minDn, " . (length($canUp) - $minUp + $dif) . "\n";
	print STDERR "\# REF: "  . bedToStr(\@{$ref}) . "\n";
	print STDERR "\# CAN: "  . bedToStr(\@{$can}) . "\n";

	my $canSeq = substr($canUp, length($canUp) - $minUp + $dif, $minUp - $dif) . ":" . $canVar . ":" .
	    substr($canDn, 0, $minDn + $dif);
	my $refSeq = substr($refUp, length($refUp) - $minUp, $minUp) . ":" . $refVar . ":" .
	    substr($refDn, 0, $minDn);
	print STDERR "\# CAN: $canSeq\n";
	print STDERR "\# REF: $refSeq\n";
	exit(1);
    }

    ## Ensure coordiate difference update did not make downstream
    ##  buffer shorter than one sequence
    if (length($canDn) < $minDn + $dif) {
	my $dif2 = length($canDn) - ($minDn + $dif);
	$minDn += $dif2;
    }

    my $canSeq  = "";
    my $canSeqR = "";
    if (length($canUp) - $minUp + $dif < length($canUp)) {
	$canSeq  = substr($canUp, length($canUp) - $minUp + $dif, $minUp - $dif) . $canVar . 
	    substr($canDn, 0, $minDn + $dif);
	$canSeqR = substr($canUp, length($canUp) - $minUp + $dif, $minUp - $dif) . ":" . $canVar . ":" .
	    substr($canDn, 0, $minDn + $dif);
    } else {
	$canSeq  = $canVar . substr($canDn, 0, $minDn + $dif);
	$canSeqR = ":" . $canVar . ":" . substr($canDn, 0, $minDn + $dif);
    }

    my $refSeq  = substr($refUp, length($refUp) - $minUp, $minUp) . $refVar . 
	substr($refDn, 0, $minDn);
    my $refSeqR = substr($refUp, length($refUp) - $minUp, $minUp) . ":" . $refVar . ":" .
	substr($refDn, 0, $minDn);

    $canSeq =~ s/-//gi;
    $refSeq =~ s/-//gi;

    $canSeqR =~ s/-//gi;
    $refSeqR =~ s/-//gi;

    @ret = (TRUE, $canSeqR, $refSeqR);

    return(@ret) if ($canSeq eq $refSeq);

    @ret = (FALSE, $canSeqR, $refSeqR);
    return(@ret);
}

#------------------------------------------------------------------------------
sub loadVcf($$$) {
    my ($file, $dat, $header) = @_;

    mssg("\# Loading VCF file: $file.\n");

    open(VCF, "gunzip -fc $file |")
	or die("Failed to open $file for reading: $!");
    while(<VCF>) {
	chomp();

	my $line = $_;
	if ($line =~ m/^\#/) {
	    push @{$header}, $line;
	} else {
	    my @fields = split(/\t/, $line);
	    my $len = abs(getSvLen($fields[VCF_INFO]));
	    next if ($len < $options{'vcfLenThresh'});
	    $$dat{$fields[VCF_CHROM]}{$fields[VCF_POS]} = [ @fields ];
	    #push @{$dat}, [ @fields ];
	}
    }
    close(VCF);
    mssg("\# Done loading " . scalar(@{$dat}) . " VCF entries.\n\n");
}

#------------------------------------------------------------------------------
sub writeVcfRecord($$) {
    my ($fh, $dat) = @_;

    my $str = "";
    foreach my $val (@{$dat}) {
	if (defined($val)) {
	    $str .= $val;
	}
	$str .= "\t";
    }
    $str =~ s/\s+$//;
    print $fh $str . "\n";
}

#------------------------------------------------------------------------------
sub isValidType($$) {
    my ($type1, $type2) = @_;

    ## This is out of date now and should not be used given proper BKPTs now

    return(TRUE) if ($type1 =~ m/BND/ || $type2 =~ m/BND/);

    ## Remove length identifier for insertions
    $type1 =~ s/^(INS).*/$1/;
    $type2 =~ s/^(INS).*/$1/;

    return(TRUE) if ($type1 eq $type2);
    return(FALSE);
}

#------------------------------------------------------------------------------
sub getRepOverlap($) {
    my ($dat) = @_;

    my $type1 = $$dat[INT_TYPE1];
    my $type2 = $$dat[INT_TYPE2];

    ## Can't measure break point lengths
    return(TRUE) if ($type1 =~ m/BND/ || $type2 =~ m/BND/);

    my ($s1, $e1, $s2, $e2) = ($$dat[INT_START1], $$dat[INT_END1],
			       $$dat[INT_START2], $$dat[INT_END2]);

    my $len1 = abs($e1 - $s1);
    my $len2 = abs($e2 - $s2);

    my $max = $len1;
    $max = $len2 if ($len2 > $len1);

    if ($$dat[INT_TYPE1] =~ m/INS:([0-9]+)$/) { 
	$len1 = $1;
	if ($$dat[INT_TYPE2] =~ m/INS:([0-9]+)$/) {
	    $len2 = $1;

	    if (! defined($len1) || ! defined($len2)) {
		print STDERR "[ERROR]: Insertion length extraction failure: '$len1', '$len2'\n";
		for (my $ii = 0; $ii < @{$dat}; $ii++) {
		    print STDERR "\t$ii\t'$$dat[$ii]'\n";
		}
		exit(1);
	    }

	    $max = $len1;
	    $max = $len2 if ($len2 > $len1);

	    return(($max - abs($len2 - $len1)) / $max);
	}
    }

    ## Otherwise find the largest start and lowest end to calculate area
    my $sM = $s1;
    if ($s2 > $sM) { $sM = $s2; }
    my $eM = $e1;
    if ($e2 < $eM) { $eM = $e2; }

    my $area = $eM - $sM + 1;

    my $percent = $area / $max;
    #if ($area / $max == 0.693418940609952) {
    if ($percent =~ m/0.69341894/) {
	print STDERR "Values: $s1, $e1, $s2, $e2, LEN: $len1, $len2, AREA: $area, $max, SM: $sM, $eM\n";
    }

    return($area / $max);
}

#------------------------------------------------------------------------------
sub getCall($) {
    my ($dat) = @_;

    ##
    ## 1 = >90\% reciprocal overlap and same type
    ## 2 = <90\% reciprocal overlap and same type
    ## 3 = break point call type
    ## 4 = mismatched type
    ## 4 = no match
    ##

    #print STDERR "C: @{$dat}\n";

    ## First take care of the break poitn cases
    return(PASS_CALL) if ($$dat[INT_TYPE1] =~ m/BND/ && $$dat[INT_TYPE2] =~ m/BND/);
    return(BKPT_CALL) if ($$dat[INT_TYPE1] =~ m/BND/ || $$dat[INT_TYPE2] =~ m/BND/);

    ## Quick validation for Grouper format updates
    if (! defined($$dat[INT_TYPE1]) || ! defined($$dat[INT_TYPE2])) {
	print STDERR "[ERROR]: Failed before isValidType: '$$dat[INT_TYPE1]', or '$$dat[INT_TYPE2]'\n";
	foreach my $val (@{$dat}) {
	    print STDERR "Data value: $val\n";
	}
	die("Failed some new VCF update spec....");
    }

    ## First check if type is correct
    my $validType = isValidType($$dat[INT_TYPE1], $$dat[INT_TYPE2]);

    return(MISTYPE_CALL) if (! $validType);

    my $overlap = getRepOverlap(\@{$dat});
    #print STDERR "$overlap >= $options{'rcOverlap'}\n";

    return(PASS_CALL) if ($overlap >= $options{'rcOverlap'});
    return(MISRP_CALL);
}

#------------------------------------------------------------------------------
sub getSvLen($) {
    my ($info) = @_;

    my $len;
    if ($info =~ m/SVLEN_C=([-0-9]+)/) {
	$len = $1;
    } else {
	die("Failed to parse SV Length from: $info!");
    }
    return($len);
}

#------------------------------------------------------------------------------
sub getSvType($) {
    my ($info) = @_;

    my $len;
    if ($info =~ m/SVTYPE_C=([A-Z:]+)/) {
	$len = $1;
    } else {
	die("Failed to parse SV Length from: $info!");
    }
    return($len);
}

#------------------------------------------------------------------------------
sub writeStats($$$) {
    my ($statType, $statFile, $d) = @_;

    my %dat = %{$d};

	print ("\# Stat results for $statType: $statFile\n");	
	foreach my $svType (@svTypes) {
	    print ("$svType:\nSV_Size");
	    foreach my $svCall (@svCalls) {
		print ("\t$svCall");
	    }
	    print ("\n");

	    #foreach my $svBin (@ranges) {
	    for (my $ii = 0; $ii < @ranges; $ii++) {
		my $sum = 0;
		foreach my $svCall (@svCalls) {
		    if (defined($dat{$svType}{$ii}{$svCall})) {
			$sum += $dat{$svType}{$ii}{$svCall}
		    }
		}
		my $svBin = $ranges[$ii];
		$svBin =~ s/,/ to /;
		print ("$svBin ($sum)");
		if ($sum == 0) {
		    foreach my $svCall (@svCalls) {
			print ("\t0 (0)");
		    }
		} else {
		    foreach my $svCall (@svCalls) {
			my $val = "0";
			if (defined($dat{$svType}{$ii}{$svCall})) {
			    $val = $dat{$svType}{$ii}{$svCall} / $sum;
			    $val = sprintf("%.2f", $val);
			    $val = $val * 100;
			    print ("\t$val% ($dat{$svType}{$ii}{$svCall})"); 
			} else {
			    print ("\t$val% (0)");
			}
		    }
		}
		print ("\n");
	    }
	    print ("\n");
	}
	print ("\n");
}

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#
# Program Execute
#
#------------------------------------------------------------------------------

$params{'start-time'} = time();

my %refDat = ();
my %refRem = ();

processArgs();

## First load all reference targets
my @refSources = ();

## Process each candidate target and produce summaries
my $intIdx = 0;
my $intCanIdx = 0;
foreach my $canVcfFile (@{$options{'can-files'}}) {
    ## Init candidate data structures
    my %canDat    = ();
    my @canHeader = ();
    my $canSrc    = loadVcf($canVcfFile, \@canDat, \@canHeader);

    ##
    ## Now we will process each intersection against each reference
    ##  We will need to keep track of the order of the references
    ##  to update VCF header of the candidate file
    ##
    my $intRefIdx = 1;
    foreach my $refVcfFile (@{$options{'ref-files'}}) {
	open(REFIN, "<$refVcfFile")
	    or die("Failed to open $refVcfFile for reading: $!");
	mssg("\# Processing refFile: $refVcfFile\n");
	my $refIdx = 0;
	while(<REFIN>) {
	    chomp();
	    my $line = $_;
	    
	    ## Handle VCF header
	    if ($line =~ m/^\#\#/) {
		print REFOUT "$line\n";
		next;
	    } elsif ($line =~ m/^\#/) {
		print REFOUT "$refHeader";
		print REFOUT "$line\n";
		next;
	    }

	    my @refDat = split(/\t/, $line);

	    my $bestCall = WORST_CALL;
	    my @canIndicies = ();
	    if (defined($refCalls[$refIdx])) {
		my @tmp = sort {$a <=> $b} keys %{$refCalls[$refIdx]};
		foreach my $key (@tmp) {
		    $bestCall = $key if ($key < $bestCall);
		}
	    }

	    ## Record stats
	    my $svLen  = abs(getSvLen($refDat[VCF_INFO]));
	    my $svType = getSvType($refDat[VCF_INFO]);
	    my $svBin  = $map{$svLen};
	    if ( ! defined($svBin)) {
		$svBin = $lastBin;
	    }
	    if (! defined($refCallStats{$svType}{$svBin}{$bestCall})) {
		$refCallStats{$svType}{$svBin}{$bestCall} = 1;
	    } else {
		$refCallStats{$svType}{$svBin}{$bestCall}++;
	    }

	    my $matchSrc   = "";
	    my $matchType  = "";
	    my $matchIdx   = "";
	    my @bestIds = ();
	    if ($bestCall != WORST_CALL) {
		foreach my $idx (sort {$a <=> $b} @{$refCalls[$refIdx]{$bestCall}}) {
		    $matchSrc  .= "1,";
		    $matchType .= "$bestCall,";
		    $matchIdx  .= "$idx,";
		    push @bestIds, $idx;
		}
	    } else {
		$matchSrc  = 1;
		$matchType = $bestCall;
		$matchIdx  = -1;
	    }
	    $matchSrc  =~ s/,+$//;
	    $matchType =~ s/,+$//;
	    $matchIdx  =~ s/,+$//;

	    if (! defined($refDat[VCF_INFO])) {
		die("\# No info field line: '$line'");
	    }

	    if (FALSE && $refDat[1] == 42773985) {
		foreach my $id (@bestIds) {
		    print "$refIdx MATCHSRC=$matchSrc;MATCHTYPE=$matchType;MATCHIDX=$id\n";
		    print "REF: $line\n";
		    print "CAN: @{$canDat[$id]}\n";
		}
	    }

	    ## Update Info
	    $refDat[VCF_INFO] =~ s/;+$//;
	    $refDat[VCF_INFO] .= ";MATCHSRC=$matchSrc";
	    $refDat[VCF_INFO] .= ";MATCHTYPE=$matchType";
	    $refDat[VCF_INFO] .= ";MATCHIDX=$matchIdx";

	    writeVcfRecord(\*REFOUT, \@refDat);
	    $refIdx++;
	}
	close(REFOUT);

	$intRefIdx++;
	writeStats("REF", $refVcfFile, \%refCallStats);
    }

    ## Open file handles
    my $canOutFile = "$options{'out-dir'}/can-$intCanIdx.vcf";
    open(CANOUT, ">$canOutFile")
	or die("Failed to open $canOutFile for writing: $!");
    open(CANIN, "<$canVcfFile")
	or die("Failed to open $canVcfFile for reading: $!");
    mssg("\# Processing canFile: $canVcfFile\n");
    my $canIdx = 0;
    my %canCallStats = ();
    while(<CANIN>) {
	chomp();
	my $line = $_;
	    
	## Handle VCF header
	if ($line =~ m/^\#\#/) {
	    print CANOUT "$line\n";
	    next;
	} elsif ($line =~ m/^\#/) {
	    ## Simply using extra fields from Reference headers
	    #print CANOUT "$canHeader";
	    print CANOUT "$refHeader";
	    print CANOUT "$line\n";
	    next;
	}
	
	my @canDat = split(/\t/, $line);
	
	my $bestCall = WORST_CALL;
	my $bestIdx;
	my @canIndicies = ();
	if (defined($canCalls[$intCanIdx][$canIdx])) {
	    my @tmp = sort {$a <=> $b} keys %{$canCalls[$intCanIdx][$canIdx]};
	    foreach my $key (@tmp) {
		$bestCall = $key if ($key < $bestCall);

		#if ($$key[1] < $bestCall) {
		#$bestCall = $$key[1];
		#   $bestIdx  = $$key[0];
		#}
	    }
	}

	## Record stats
	my $svLen  = abs(getSvLen($canDat[VCF_INFO]));
	my $svType = getSvType($canDat[VCF_INFO]);
	my $svBin  = $map{$svLen};
	if ( ! defined($svBin)) {
	    $svBin = $lastBin;
	}
	if (! defined($canCallStats{$svType}{$svBin}{$bestCall})) {
	    $canCallStats{$svType}{$svBin}{$bestCall} = 1;
	} else {
	    $canCallStats{$svType}{$svBin}{$bestCall}++;
	}

	my $matchSrc   = "";
	my $matchType  = "";
	my $matchIdx   = "";
	my @bestIds = ();
	## FIXME: PICK UP HERE!!!
	if ($bestCall != WORST_CALL) {
	    foreach my $idx (sort {$a <=> $b} @{$canCalls[$intCanIdx][$canIdx]{$bestCall}}) {
		my $bestIdx = $$idx[0];
		$matchSrc  .= "$bestIdx,";
		$matchType .= "$bestCall,";
		$matchIdx  .= "$$idx[1],";
		push @bestIds, $idx;
	    }
	} else {
	    $matchSrc  = 1;
	    $matchType = $bestCall;
	    $matchIdx  = -1;
	}
	$matchSrc  =~ s/,+$//;
	$matchType =~ s/,+$//;
	$matchIdx  =~ s/,+$//;
	
	if (! defined($canDat[VCF_INFO])) {
	    die("\# No info field line: '$line'");
	}
	
	if (FALSE && $canDat[1] == 42773985) {
	    foreach my $id (@bestIds) {
		print "$canIdx MATCHSRC=$matchSrc;MATCHTYPE=$matchType;MATCHIDX=$id\n";
		print "CAN: $line\n";
		print "CAN: @{$canDat[$id]}\n";
	    }
	}
	
	## Update Info
	$canDat[VCF_INFO] =~ s/;+$//;
	$canDat[VCF_INFO] .= ";MATCHSRC=$matchSrc";
	$canDat[VCF_INFO] .= ";MATCHTYPE=$matchType";
	$canDat[VCF_INFO] .= ";MATCHIDX=$matchIdx";
	
	writeVcfRecord(\*CANOUT, \@canDat);
	$canIdx++;
    }
    close(CANIN);
    close(CANOUT);
    mssg("\# Processed $canIdx candidates from VCF file!\n");

    writeStats("CAN", $canVcfFile, \%canCallStats);
}


## Write summary 
$params{'end-time'} = time();
$params{'run-time'} = ($params{'start-time'} - $params{'end-time'}) / 60;

if ($params{'run-time'} < 60) {
    $params{'run-time'} = sprintf("%.2f", $params{'run-time'}) . " minutes";
} else {
    $params{'run-time'} /= 60;
    $params{'run-time'} = sprintf("%.2f", $params{'run-time'}) . " hours";
}

#------------------------------------------------------------------------------

## End of file
