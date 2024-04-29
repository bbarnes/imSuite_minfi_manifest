#! /usr/bin/perl -w

use strict;
use Getopt::Long;

use constant FALSE => 0;
use constant TRUE  => 1;

my $revComp = FALSE;
my $comp    = FALSE;
my $rev     = FALSE;
my $name;
my $file;
#my $start = 48000000;
#my $len   =  2000000;
my $start;
my $len;
my $end;


my $buf = 0;
my $chr;
my $bedFile;

GetOptions(
	   "f=s"   => \$file,
	   "n=s"   => \$name,
	   "s=i"   => \$start,
	   "l=i"   => \$len,
	   "e=i"   => \$end,
	   "rc"    => \$revComp,
	   "r"     => \$rev,
	   "c"     => \$comp,
	   "b=i"   => \$buf,
	   "chr=s" => \$chr,
	   "bed=s" => \$bedFile,
	   );

if (! defined($bedFile) || ! defined($chr)) {
    print STDERR "Usage: $0 -bed bedFile -f file -chr chromosome -n name -s start [-l len or -e end] [options]\n";
    print STDERR "Usage: Missing arguments! Exiting...\n";
    exit(1);
}

#    ! defined($start) ||
#    (! defined($len) && ! defined($end)) ||
#    ! defined($name) ||
#    ! defined($chr) ||


die("Start > end: $start > $end!") if (defined($end) && $end < $start);

$file = "/illumina/scratch/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/$chr.fa";

open(IN, "<$file")
    or die "Could not open $file for reading: $!";
my $genome = "";
my $lineCnt = 0;
while(<IN>) {
    s/^\s+//;
    s/\s+$//;
    next if /^(\s)*$/;
    next if /^>.*/;
    
    $lineCnt++;
    
    my $seq = $_;
    $genome .= $seq;
}
close(IN);

print STDERR "Len=" . length($genome) . "\n";

open(BED, "<$bedFile")
    or die("Failed to open $bedFile for reading!");

while(<BED>) {
    s/^\s+//;
    s/\s+$//;
    next if /^\s*$/;

    my $line = $_;
    my @fields = split(/\t/, $line);
    my $curChr = $fields[0];
    my $start  = $fields[1];
    my $end    = $fields[2];
    my $name   = $fields[3];

    # Normalize for genome subseq
    #$start--;

    $start -= $buf;

    if (! defined($name)) { $name = "$curChr.$start-$end"; }
    else { $name = "$name:$curChr.$start-$end"; }

    if (! defined($len)) {
	$len = $end - $start + 1 + $buf;
    }
    
    my $seq = substr($genome, $start, $len);
    
    if ($revComp) {
	$seq = reverse($seq);
	$seq =~ tr/actgACTG/tgacTGAC/;
    }
    if ($rev) {
	$seq = reverse($seq);
    }
    if ($comp) {
	$seq =~ tr/actgACTG/tgacTGAC/;
    }

    print ">$name\n$seq\n";
}
close(BED);

## End of file
