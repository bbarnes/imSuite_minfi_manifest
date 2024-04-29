#! /usr/bin/perl -w

use strict;
use Getopt::Long;

use constant FALSE => 0;
use constant TRUE  => 1;

## Bed file formats
use constant BED_CHR   => 0;
use constant BED_START => 1;
use constant BED_END   => 2;
use constant BED_ID    => 3;
use constant BED_IDX   => 4;
use constant BED_SIDE  => 5;
use constant BED_INFO  => 6;

my @mapDat  = ();
my @infoDat = ();
my @results = ();

my $numElements;
my $sampleName = "NA";

GetOptions("n=i"  => \$numElements,
	   "s=s"  => \$sampleName,
	   );

if (! defined($numElements) || ! defined($sampleName)) {
    print STDERR "[Usage]: $0 -n numberElements -s sampleName [options]\n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n";
    exit(1);
}

my %lineHash1 = ();
my %lineHash2 = ();
my $lineCnt = 0;
while(<STDIN>) {
    s/^\s+//;
    s/\s+$//;
    next if /^\#/;

    my $line = $_;
    my @fields = split(/\t/, $line);

    #print "LINE: $line\n";

    my @idxes = split(/,/, $fields[BED_IDX]);
    my @sides = split(/,/, $fields[BED_SIDE]);
    my @infos = split(/@@/, $fields[BED_INFO]);

    ## Loop over samples
    for (my $smpIdx = 0; $smpIdx < @idxes; $smpIdx++) {
	my $idx  = $idxes[$smpIdx];
	my $side = $sides[$smpIdx];
	my $info = $infos[$smpIdx];

	## Skip if variant is not found in sample
	next if ($idx eq "-");

	## 
	## Sanity check: Skip for now since they can now how have multiple values
	##
	die("[ERROR]: [Line: $lineCnt]: Failed to match info ($smpIdx), line:\n" .
	    "$line\n\t'$infoDat[$smpIdx][$idx]'\n\t'$info'\n\n" .
	    "\tHASH1\t" .   $lineHash1{"$smpIdx.$idx.$side"} . 
	    "\n\tHASH2\t" . $lineHash2{"$smpIdx.$idx"} . "\n" )
	    if (defined($infoDat[$smpIdx][$idx]) && $infoDat[$smpIdx][$idx] ne $info && 0);

	#die("[ERROR]: [Line: $lineCnt]: Mapping already exists: $smpIdx, $idx, $side: @{$mapDat[$smpIdx][$idx][$side]}!\nOn line: '$line'\n" .
	#   "$infos[$smpIdx]\n$infoDat[$smpIdx][$idx]\n" . "\tHASH1\t" . $lineHash1{"$smpIdx.$idx.$side"} . "\n\tHASH2\t" . $lineHash2{"$smpIdx.$idx"} . "\n")
	#  if (defined($mapDat[$smpIdx][$idx][$side]));


	## Assign infom fields (we should only have to do this here)
	$infoDat[$smpIdx][$idx] = $infos[$smpIdx];

	## Likely junk to remove
	#print STDERR "$smpIdx\t$idx\t$infoDat[$smpIdx][$idx]\n";
	$lineHash1{"$smpIdx.$idx.$side"} = "$lineCnt:$line";
	$lineHash2{"$smpIdx.$idx"} = $line;

	## Sanity check
	die("[ERROR]: [Line: $lineCnt, smpIdx]: Failed to match info ($smpIdx), line:\n" .
	    "$line\n\t'$infoDat[$smpIdx][$idx]'\n\t'$info'\n\n" .
	    "\tHASH1\t" .   $lineHash1{"$smpIdx.$idx.$side"} . 
	    "\n\tHASH2\t" . $lineHash2{"$smpIdx.$idx"} . "\n" )
	    if ($smpIdx eq "-");

	die("[ERROR]: [Line: $lineCnt, idx]: Failed to match info ($smpIdx), line:\n" .
	    "$line\n\t'$infoDat[$smpIdx][$idx]'\n\t'$info'\n\n" .
	    "\tHASH1\t" .   $lineHash1{"$smpIdx.$idx.$side"} . 
	    "\n\tHASH2\t" . $lineHash2{"$smpIdx.$idx"} . "\n" )
	    if ($idx eq "-");

	die("[ERROR]: [Line: $lineCnt, side]: Failed to match info ($smpIdx), line:\n" .
	    "$line\n\t'$infoDat[$smpIdx][$idx]'\n\t'$info'\n\n" .
	    "\tHASH1\t" .   $lineHash1{"$smpIdx.$idx.$side"} . 
	    "\n\tHASH2\t" . $lineHash2{"$smpIdx.$idx"} . "\n" . "TRIPLET: $smpIdx, $idx, $side\n")
	    if ($side eq "-");

	## Allow multiple mappings
	push @{$mapDat[$smpIdx][$idx][$side]}, [ $fields[BED_IDX], $fields[BED_SIDE] ];
	#$mapDat[$smpIdx][$idx][$side] = [ $fields[BED_IDX], $fields[BED_SIDE] ];
    }
    #print "-" x 100 . "\n";
    $lineCnt++;
}

my %summary = ();
my $blankStr = "-," x $numElements;
$blankStr =~ s/,+$//;
for (my $smpIdx = 0; $smpIdx < @mapDat; $smpIdx++) {
    if (! defined($mapDat[$smpIdx])) {
	#print STDERR "[Warning]: Failed to find sample index for $smpIdx!\n";
	next;
    }
    for (my $idx = 0; $idx < @{$mapDat[$smpIdx]}; $idx++) {
	my @bestIndexA = split(/,/, $blankStr);
	my @bestIndexB = split(/,/, $blankStr);
	my @bestSidesA = split(/,/, $blankStr);
	my @bestSidesB = split(/,/, $blankStr);
	my $bestMatCnt = -1;
	my $bestMatStr = "";

	##
	## The first two cases below describe an event with only on side matching.
	##  For now we will report these as a not hit. However, we could come up with
	##  some more informative reporting later. 
	##
	if (defined($mapDat[$smpIdx][$idx][0]) && ! defined($mapDat[$smpIdx][$idx][1])) {
	    #print STDERR "[Warning]: A: No matches for either side (sample, index, n-elements): " .
	    #"$smpIdx, $idx, $numElements\n";
	    #print STDERR "[Warning]: @{$mapDat[$smpIdx][$idx]}\n";
	    next;
	    @bestIndexA = split(/,/, $mapDat[$smpIdx][$idx][0][0][0]);
	    $bestMatCnt = 0;
	    for (my $jj = 0; $jj < @bestIndexA; $jj++) {
		$bestMatStr .= "$jj," if ($bestIndexA[$jj] ne "-");
	    }
	} elsif (! defined($mapDat[$smpIdx][$idx][0]) && defined($mapDat[$smpIdx][$idx][1])) {
	    #print STDERR "[Warning]: B: No matches for either side (sample, index, n-elements): " .
	    #"$smpIdx, $idx, $numElements\n";
	    #print STDERR "[Warning]: @{$mapDat[$smpIdx][$idx]}\n";

	    next;
	    @bestIndexA = split(/,/, $mapDat[$smpIdx][$idx][1][0][0]);
	    $bestMatCnt = 0;
	    for (my $jj = 0; $jj < @bestIndexA; $jj++) {
		$bestMatStr .= "$jj," if ($bestIndexA[$jj] ne "-");
	    }
	} elsif (! defined($mapDat[$smpIdx][$idx][0]) && ! defined($mapDat[$smpIdx][$idx][1])) {
	    #print STDERR "[Warning]: C: No matches for either side (sample, index, n-elements): " .
	    #"$smpIdx, $idx, $numElements\n";
	    #print STDERR "[Warning]: @{$mapDat[$smpIdx][$idx]}\n";
	    next;
	} else {
	    for (my $sd1 = 0; $sd1 < @{$mapDat[$smpIdx][$idx][0]}; $sd1++) {
		my @indexA = split(/,/, $mapDat[$smpIdx][$idx][0][$sd1][0]);
		my @sidesA = split(/,/, $mapDat[$smpIdx][$idx][0][$sd1][1]);
		
		for (my $sd2 = 0; $sd2 < @{$mapDat[$smpIdx][$idx][1]}; $sd2++) {
		    my @indexB = split(/,/, $mapDat[$smpIdx][$idx][1][$sd2][0]);
		    my @sidesB = split(/,/, $mapDat[$smpIdx][$idx][1][$sd2][1]);
		    
		    ## Investigate match and decide how good it is
		    my $matCnt = 0;
		    my $matStr = "";
		    for (my $jj = 0; $jj < @indexA; $jj++) {
			#my $tmpFound = "FALSE";
			if ($indexA[$jj] eq $indexB[$jj] && $sidesA[$jj] ne $sidesB[$jj]) {
			    $matCnt++;
			    $matStr .= "$jj,";
			    #print STDERR "HERE: matStr: $matStr\n";
			    #$tmpFound = "TRUE";
			}
			#print STDERR "$smpIdx, $idx, $sd1, $sd2, $jj: $tmpFound\n";
			#print STDERR "\t@indexA, @sidesA\n";
			#print STDERR "\t@indexB, @sidesB\n";
			#print STDERR "\t$indexA[$jj] eq $indexB[$jj] && $sidesA[$jj] ne $sidesB[$jj]\n";
			#print STDERR "\tLINE1A: " . $lineHash1{"$smpIdx.$idx.$sidesA[$jj]"} . "\n";
			#print STDERR "\tLINE1B: " . $lineHash1{"$smpIdx.$idx.$sidesB[$jj]"} . "\n";
		    }
		    if ($matCnt > $bestMatCnt) {
			@bestIndexA = @indexA;
			@bestIndexB = @indexB;
			@bestSidesA = @sidesA;
			@bestSidesB = @sidesB;
			$bestMatCnt = $matCnt;
			$bestMatStr = $matStr;
		    }
		}
	    }
        }
	$bestMatStr =~ s/,+$//;

	my $outStr = "";
	my $checkStr = "";
	my @tmpIndex = split(/,/, $bestMatStr);
	my %tmpIdxHash = ();
	foreach my $val (@tmpIndex) {
	    $tmpIdxHash{$val} = 1;
	}

	for (my $jj = 0; $jj < $numElements; $jj++) {
	    #if ($bestIndexA[$jj] eq "-") {
	    if (! defined($tmpIdxHash{$jj})) {
		$outStr .= "-\t";
	    } else {
		$outStr .= "$infoDat[$jj][$bestIndexA[$jj]]\t";
		$checkStr .= "$jj,";
	    }
	}
	$checkStr =~ s/,+$//;
	if ($checkStr ne $bestMatStr) {
	    die("[ERROR]: Check string != best string: $checkStr != $bestMatStr\n");
	}

	#$outStr = "HIT\t$sampleName\t$checkStr\t" . $outStr;
	$outStr = "HIT\t$sampleName\t$bestMatStr\t" . $outStr;
	print "$outStr\n";
	#print STDERR "\n\n";
    }
}

#print STDERR "\# Done! Size, mapDat: " . @mapDat . "\n";

##
## We currently don't believe these stats
##

#foreach my $key (sort keys %summary) {
#    print "$key\t$summary{$key}\n";
#}

## End of file
