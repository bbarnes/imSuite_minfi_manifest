#! /usr/bin/perl

use strict;
use Getopt::Long;

use constant SEQ_ID                              => 0;
use constant FORWARD_SEQUENCE                    => 1;
use constant GENOME_BUILD                        => 2;
use constant CHROMOSOME                          => 3;
use constant COORDINATE                          => 4;
use constant DESIGN_STATE                        => 5;
use constant SEQ_LENGTH                          => 6;
use constant FORWARD_CPG_COORD                   => 7;
use constant TB_STRAND                           => 8;
use constant TOP_SEQUENCE                        => 9;
use constant TOP_CPG_COORD                       => 10;
use constant PROBE_TYPE                          => 11;
use constant PROBESET_ID                         => 12;
use constant PROBESET_SCORE                      => 13;
use constant METHYL_PROBE_ID                     => 14;
use constant METHYL_PROBE_SEQUENCE               => 15;
use constant METHYL_PROBE_LENGTH                 => 16;
use constant METHYL_START_COORD                  => 17;
use constant METHYL_END_COORD                    => 18;
use constant METHYL_PROBE_COVERED_TOP_SEQUENCE   => 19;
use constant METHYL_ALLELE_FR_STRAND             => 20;
use constant METHYL_ALLELE_TB_STRAND             => 21;
use constant METHYL_ALLELE_CO_STRAND             => 22;
use constant METHYL_ALLELE_TYPE                  => 23;
use constant METHYL_FINAL_SCORE                  => 24;
use constant METHYL_TM                           => 25;
use constant METHYL_TM_SCORE                     => 26;
use constant METHYL_GC_PERCENT                   => 27;
use constant METHYL_GC_SCORE                     => 28;
use constant METHYL_13MER_COUNT                  => 29;
use constant METHYL_13MER_SCORE                  => 30;
use constant METHYL_ADDRESS_COUNT                => 31;
use constant METHYL_ADDRESS_SCORE                => 32;
use constant METHYL_SELF_COMPLEMENTARITY         => 33;
use constant METHYL_SELF_COMPLEMENTARITY_SCORE   => 34;
use constant METHYL_MONO_RUN                     => 35;
use constant METHYL_MONO_RUN_SCORE               => 36;
use constant METHYL_ECTOPIC_COUNT                => 37;
use constant METHYL_ECTOPIC_SCORE                => 38;
use constant METHYL_UNDERLYING_CPG_COUNT         => 39;
use constant METHYL_UNDERLYING_CPG_MIN_DIST      => 40;
use constant METHYL_UNDERLYING_CPG_SCORE         => 41;
use constant METHYL_IN_CPG_ISLAND_RELAXED        => 42;
use constant METHYL_CPG_ISLAND_SCORE             => 43;
use constant METHYL_NEXT_BASE                    => 44;
use constant METHYL_NEXT_BASE_SCORE              => 45;
use constant UNMETHYL_PROBE_ID                   => 46;
use constant UNMETHYL_PROBE_SEQUENCE             => 47;
use constant UNMETHYL_PROBE_LENGTH               => 48;
use constant UNMETHYL_START_COORD                => 49;
use constant UNMETHYL_END_COORD                  => 50;
use constant UNMETHYL_PROBE_COVERED_TOP_SEQUENCE => 51;
use constant UNMETHYL_ALLELE_FR_STRAND           => 52;
use constant UNMETHYL_ALLELE_TB_STRAND           => 53;
use constant UNMETHYL_ALLELE_CO_STRAND           => 54;
use constant UNMETHYL_ALLELE_TYPE                => 55;
use constant UNMETHYL_FINAL_SCORE                => 56;
use constant UNMETHYL_TM                         => 57;
use constant UNMETHYL_TM_SCORE                   => 58;
use constant UNMETHYL_GC_PERCENT                 => 59;
use constant UNMETHYL_GC_SCORE                   => 60;
use constant UNMETHYL_13MER_COUNT                => 61;
use constant UNMETHYL_13MER_SCORE                => 62;
use constant UNMETHYL_ADDRESS_COUNT              => 63;
use constant UNMETHYL_ADDRESS_SCORE              => 64;
use constant UNMETHYL_SELF_COMPLEMENTARITY       => 65;
use constant UNMETHYL_SELF_COMPLEMENTARITY_SCORE => 66;
use constant UNMETHYL_MONO_RUN                   => 67;
use constant UNMETHYL_MONO_RUN_SCORE             => 68;
use constant UNMETHYL_ECTOPIC_COUNT              => 69;
use constant UNMETHYL_ECTOPIC_SCORE              => 70;
use constant UNMETHYL_UNDERLYING_CPG_COUNT       => 71;
use constant UNMETHYL_UNDERLYING_CPG_MIN_DIST    => 72;
use constant UNMETHYL_UNDERLYING_CPG_SCORE       => 73;
use constant UNMETHYL_IN_CPG_ISLAND_RELAXED      => 74;
use constant UNMETHYL_CPG_ISLAND_SCORE           => 75;
use constant UNMETHYL_NEXT_BASE                  => 76;
use constant UNMETHYL_NEXT_BASE_SCORE            => 77;

my $designFile;
my $listFile;

my $idxName = 500;

my $max;
my $zeroBased = 0;

GetOptions("d=s"   => \$designFile,
	   "l=s"   => \$listFile,
	   "idx=i" => \$idxName,
	   "max=i" => \$max,
	   );

if (! defined($designFile) ) {
    print STDERR "[Usage]: $0 -d designFile [ -l listFile ] -idx idxName [options]\n";
    print STDERR "[Usage]: Missing argument! Exiting...\n";
    exit(1);
}

my %targets = ();
if (defined($listFile)) {
    open(IN, "<$listFile")
	or die("Failed to open $listFile for reading: $!");
    while(<IN>) {
	s/^\s+//;
	s/\s+$//;

	my @fields = split(/\t/, $_);
	$targets{$fields[0]} = 1;
    }
    close(IN);
    print STDERR "\# Loaded " . scalar(keys %targets) . " total target CGs!\n";
}

print STDERR "[Input]: designFile=$designFile!\n";

if ($designFile =~ m/\.gz/) {
    open(IN, "gzip -dc $designFile | ")
	or die("Failed to open gzip -dc $designFile for reading: $!");
} else {
    open(IN, "<$designFile")
	or die "Could not open $designFile for reading; $!";
}

my $totCnt = 0;
my $cnt = 0;
while(<IN>) {
    s/^\s+//;
    s/\s+$//;
    next if /^(\s)*$/;
    next if /^\#.*$/;
    next if /^Seq_ID.*$/;

    $totCnt++;

    my $line = $_;
    my @data = split(/\t/, $line);

    my $cgn = $data[SEQ_ID];
    my $seq = $data[FORWARD_SEQUENCE];
    my $bld = $data[GENOME_BUILD];
    my $chr = $data[CHROMOSOME];
    my $pos = $data[COORDINATE];

    my $fr  = $data[METHYL_ALLELE_FR_STRAND];
    my $co  = $data[METHYL_ALLELE_CO_STRAND];
    my $tb  = substr($data[METHYL_ALLELE_TB_STRAND],0,1);

    next if ($fr ne "F");
    next if ($co ne "C");
    next if ($tb ne "T");

    # Could encode srd (TB,CO,FR) as bit...
    my $srd_bit=-1;
    if (0) {
	$srd_bit=-1;
	
    } elsif ($tb eq "T" && $co eq "C" && $fr eq "F") {
	$srd_bit=0;
    } elsif ($tb eq "T" && $co eq "O" && $fr eq "F") {
	$srd_bit=1;
    } elsif ($tb eq "B" && $co eq "C" && $fr eq "F") {
	$srd_bit=2;
    } elsif ($tb eq "B" && $co eq "O" && $fr eq "F") {
	$srd_bit=3;

    } elsif ($tb eq "T" && $co eq "C" && $fr eq "R") {
	$srd_bit=4;
    } elsif ($tb eq "T" && $co eq "O" && $fr eq "R") {
	$srd_bit=5;
    } elsif ($tb eq "B" && $co eq "C" && $fr eq "R") {
	$srd_bit=6;
    } elsif ($tb eq "B" && $co eq "O" && $fr eq "R") {
	$srd_bit=7;

    } else {
	$srd_bit-1;
    }
    
    # Convert cg to integer
    my $cgn_int=$cgn;
    $cgn_int =~ s/^cg//;
    $cgn_int =~ s/^0+//;

    # Clean chromosome
    my $chr_int = $chr;
    $chr_int =~ s/^chr//;
    $chr_int = uc($chr_int);
    if ($chr_int =~ m/^X/) {
	$chr_int = 23;
    } elsif ($chr_int =~ m/^Y/) {
	$chr_int = 24;
    } elsif ($chr_int =~ m/^M/) {
	$chr_int = 25;
    }

    if ($chr_int !~ m/^[0-9]+$/) {
	print STDERR "[ERROR]: Failed to convert chr to int: $chr_int!\n";
	exit(1);
    }

    # Strip fwd_seq of brackets
    my $fwd_seq = $seq;
    $fwd_seq =~ s/\[//;
    $fwd_seq =~ s/\]//;

    print join("\t", $cgn_int,$srd_bit,$chr_int,$pos,$fwd_seq)."\n";
    $cnt++;

    last if (defined($max) && $cnt>=$max);
}
close(IN);

print STDERR "\# Proccessed $cnt probes out of $totCnt\n";

## End of file
