#! /usr/bin/perl -w

use strict;
use Getopt::Long;

use constant GENEA_1   => 1;
use constant CHRA_1    => 2;
use constant POSA_1    => 3;

use constant GENEB_1   => 4;
use constant CHRB_1    => 5;
use constant POSB_1    => 6;

use constant IS_GENIC  => 9;

my $format;

## All filters
my $skipGenic;
my $diffGenes;
my $splitThresh;
my $spanThresh;
my $distThresh;

## RNA Express Filters
my $skipFail;

## Next Bio Filters
my $nb_cluster_thresh;
my $nb_fusion_thresh;

## Isis RNA Filters
my $scoreThresh = 0;

my %allowedFormats = ("NB"    => 1,
		      "Isis"  => 1,
		      "Exp"   => 1,
		      "True"  => 1,
		      );

GetOptions("f=s"       => \$format,
	   "skipFail"  => \$skipFail,
	   "skipGenic" => \$skipGenic,
	   "diffGenes" => \$diffGenes,
	   "distThresh=i" => \$distThresh,

	   "scoreThresh=f"   => \$scoreThresh,

           "nb_cluster_thresh=f" => \$nb_cluster_thresh,
           "nb_fusion_thresh=f"  => \$nb_fusion_thresh,
           "splitThresh=i"       => \$splitThresh,
           "spanThresh=i"        => \$spanThresh,
	   );

if (! defined($format)) {
    print STDERR "[Usage]: $0 -f [NB,Isis,Exp] [options] < inputFile.txt > output.bed";
    print STDERR "[Usage]: Missing arguments! Exiting...\n";
    exit(1);
}

die("[ERROR]: Unknown format: $format")
    if (! defined($allowedFormats{$format}));

my $eventCnt = 0;
my $posCnt   = 0;
while(<STDIN>) {
    s/^\s+//;
    s/\s+$//;
    next if /^(\s)*$/;
    next if /^\#/;

    my $line = $_;
    my @fields = split(/\t/, $line);

    die("[ERROR]: Line containts seperator ':' on line: $line") if ($line =~ m/:/);

    my $squashLine = $line;
    $squashLine =~ tr/\t/:/;

    my ($chrA, $posA,
	$chrB, $posB,
	$geneA, $geneB,
	$splitReads, $spanReads,
	$nb_cluster,
	$nb_fusion,
	$score
	);
    my $dist = -1;
    if ($format eq "Isis") {
	$chrA  = $fields[CHRA_1];
	$posA  = $fields[POSA_1];
	$chrB  = $fields[CHRB_1];
	$posB  = $fields[POSB_1];
	$geneA = $fields[GENEA_1];
	$geneB = $fields[GENEB_1];
	$score = $fields[10];

	next if (defined($scoreThresh) && $score < $scoreThresh);

	$dist = abs($posA - $posB) if ($chrA eq $chrB);

	if ($fields[7] =~ m/;/) {
	    my @tmp = split(/;/, $fields[7]);
	    $splitReads = $tmp[1];
	} else {
	    $splitReads = $fields[7];
	}

	if ($fields[8] =~ m/;/) {
	    my @tmp = split(/;/, $fields[8]);
	    $spanReads  = $tmp[1];
	} else {
	    $spanReads  = $fields[8];
	}
    } elsif ($format eq "Exp") {
	$chrA  = $fields[CHRA_1];
	$posA  = $fields[POSA_1];
	$chrB  = $fields[CHRB_1];
	$posB  = $fields[POSB_1];
	$geneA = $fields[GENEA_1];
	$geneB = $fields[GENEB_1];

	$dist = abs($posA - $posB) if ($chrA eq $chrB);

	my @tmp = split(/;/, $fields[7]);
	$splitReads = $tmp[1];
	@tmp = split(/;/, $fields[8]);
	$spanReads  = $tmp[1];
    } elsif ($format eq "NB") {
	next if ($line =~ m/^cluster_score/);
	if ($fields[3] =~ m/^([^-]+)-([^-]+)$/) {
	    $chrA = $1;
	    $chrB = $2;
	} else {
	    die("[ERROR]: Unable to parse NB results; $fields[3], on line: $line!\n");
	}
	$posA = $fields[4];
	$posB = $fields[5];

	if ($fields[2] =~ m/^([^-]+)-([^-]+)$/) {
	    $geneA = $1;
	    $geneB = $2;
	} elsif ($fields[2] =~ m/^([^-]+-[^-]+)-([^-]+-[^-]+)$/) {
	    $geneA = $1;
	    $geneB = $2;
	} elsif ($fields[2] =~ m/^([^-]+-[^-]+)-([^-]+)$/) {
	    $geneA = $1;
	    $geneB = $2;
	} else {
	    die("[ERROR]: Unable to parse genes from NB data: $fields[2] on line: $line!");
	}

	$splitReads = $fields[8];
	$spanReads  = $fields[7];

	$nb_cluster = $fields[0];
	$nb_fusion  = $fields[1];

	$dist = abs($posA - $posB) if ($chrA eq $chrB);
	#print STDERR "$geneA, $geneB, $posA, $posB, $dist, $chrA, $chrB\n";
    } else {
	die("[ERROR]: Format: $format not supported!");
    }
    next if ($skipFail && $line =~ m/\tFalse$/);
    #next if ($chrA eq $chrB);
    if ($fields[IS_GENIC] =~ m/^[^;]+;[^;]+;[^;]+$/) { 
	#print STDERR "HERE: $line\n";
	if ($fields[IS_GENIC] =~ m/^[^;]+;False;[^;]+$/) {
	    #print STDERR "Is False";
	    next if ($skipGenic);
	}
    }

    if (! defined($geneA) || ! defined($geneB)) {
	die("ERROR]: Genes not defined on line: $line");
    }

    ## Skip if both genes are the same
    if (defined($diffGenes) && $diffGenes && $geneA eq $geneB) {
	#print STDERR "SAME: '$geneA' == '$geneB'\n";
	next;
    }

    ## Read filters
    if (! defined($splitReads) || $splitReads =~ m/;/) {
	my @tmp = split(/\;/, $fields[7]);
	print STDERR "TMP($format): '@tmp'\n";
	print STDERR "[ERROR]: FAILED ON THIS LINE: $line\n";
	print STDERR "\t$fields[7], $fields[8]!\n";
	print STDERR "\t$splitReads, $spanReads\n\n";
	exit(1);
    }

    next if (defined($splitThresh) && ($splitReads eq "-" || $splitReads <= $splitThresh));
    next if (defined($spanThresh) && ($spanReads eq "-" || $spanReads <= $spanThresh));
    next if (defined($nb_cluster_thresh) && $format eq "NB" && $nb_cluster <= $nb_cluster_thresh);
    next if (defined($nb_fusion_thresh) && $format eq "NB" && $nb_fusion <= $nb_fusion_thresh);
    next if (defined($distThresh) && $dist < $distThresh && $dist != -1);

    ## Write two line one for each break point. 
    print "$chrA\t$posA\t" . ($posA + 1) . "\t$posCnt\t$eventCnt\t0\t$squashLine\n";
    $posCnt++;
    print "$chrB\t$posB\t" . ($posB + 1) . "\t$posCnt\t$eventCnt\t1\t$squashLine\n";
    $posCnt++;
    $eventCnt++;
}

##
## Subroutines
##
sub loadBed {
    my ($bedFile) = @_;

    my @dat = ();
    open(BED, "<$bedFile")
	or die("Failed to open $bedFile for reading: $!");

    while(<BED>) {
	s/^\s+//;
	s/\s+$//;
	next if /^#/;

	my $line = $_;
	my @fields = split(/\t/, $line);

	die("[ERROR]: Failed to parse 6 fields from: $line! Got " . @fields) if (@fields != 7);
	push @dat, [ @fields ];
    }	
    close(BED);
    return(@dat);
}

## End of file
