#! /usr/bin/perl -w

use strict;

my %mapM = ();
my %mapD = ();
my @dNucsM = qw(Y Y B H B H N);
my @dNucsU = qw(T T K W K W D);
my @cNucs  = qw(C Y S M B H V);
my @gNucs  = qw(G R S K B D V);
my @hNucs  = qw(A T C Y W M H N);
my ($inp, $out);
for (my $ii=0; $ii <@cNucs; $ii++) {
    my $cc = $cNucs[$ii];
    my $dd = $dNucsM[$ii];
    foreach my $gg (@gNucs) {
	$inp = $cc.$gg;
	$mapM{$inp} = lc($cc).$gg;
	$mapD{$inp} = lc($dd).$gg;
	# print ("[mapM]: $inp => ".lc($cc).$gg." [mapD]: $inp => ".lc($dd).$gg."\n");

	if (0) {
	    print "MAP_M[['$inp']] <- '".lc($cc).$gg."'\n";
	} else {
	    print "MAP_D[['$inp']] <- '".lc($dd).$gg."'\n";
	}
    }
}



# End of file
