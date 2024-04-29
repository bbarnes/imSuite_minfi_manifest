#! /usr/bin/perl -w

use strict;
use Getopt::Long;

my @vals = ();
my $max;
my %val_hash = ();

my $del = "\t";
my $prefix;
my $suffix;

my $tic = "'";

GetOptions(
	   "d=s"  => \$del,
	   "p=s"  => \$prefix,
	   "s=s"  => \$suffix,
	   );


#print "\#! /usr/bin/perl -w\n\n";
#print "";

while(<STDIN>) {
    s/^\s+//;
    s/\s+$//;
    next if /^(\s)*$/;

    my $line = $_;
    $line =~ s/\"//gi;
    #my @line_array = split(/\s+/, $line);
    my @line_array = split(/$del/, $line);

    my $ii = 0;
    foreach my $val (@line_array) {
	## Clean field
	$val =~ s/^\#+//;
	$val =~ s/^\s+//;
	$val =~ s/\s+$//;

	$val =~ s/\s/_/gi;
	$val =~ s/-/_/gi;
	$val =~ s/\)//gi;
	$val =~ s/\(//gi;
	$val =~ s/[\']//gi;
	$val =~ s/\'/_/gi;
	$val =~ s/\`//gi;
	$val =~ s/\"//gi;
	$val =~ s/\./_/gi;

	$val = "IDX" if (! defined($val) || length($val) == 0);

	my $extra = 1;
	my $orig  = $val;
        while(defined($val_hash{$val})) {
	    $val = $orig . "_$extra";
	    $extra++;
	}
	if (defined($prefix)) {
	    $val = uc($prefix) . $val;
	}
	if (defined($suffix)) {
	    $val .= uc($suffix);
	}

	$val_hash{$val} = 1;
	
	push @vals, $val;
	if (! defined($max) || $max < length($val)) {
	    $max = length($val);
	}
    }
    last;
}

for (my $ii = 0; $ii < @vals; $ii++) {
    my $val = $vals[$ii];
    print "use constant ";
    print uc($val) . " " x ($max + 1 - length($val)) . "=> $ii;\n";
}
