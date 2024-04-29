#! /usr/bin/perl -w 

use strict;
use Getopt::Long;

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

@ranges = ("1,299",
          "300,10000"
          );

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

my @qscoreRange = ([0 ,14],
		   [15,29]
		   );


while(<STDIN>) {
    s/^\s+//;
    s/\s+$//;

    my $file = $_;
    if ($file =~ m/\.gz) {
	open(IN, "gzip -dc $file | ")
	    or die("Failed to uncompress $file for reading: $!");
    } else {
	open(IN, "<$file")
	    or die("Failed to open $file for reading: $!");
    }
    while(<IN>) {
	s/^\s+//;
	s/\s+$//;
	next if /^\#/;

	my $line = $_;
	my @fields = split(/\t/, $line);

	my $qscore = $fields[5];

    }
    close(IN);
}



## End of file
