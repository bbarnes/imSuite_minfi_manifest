#! /usr/bin/perl -w

use strict;

my $binDir = "/home/bbarnes/vcf-compare/scripts";

my $rc = 0.50;

my $buildDir = "/illumina/builds/platinumgenomes/Platinum_Genomes_Short_Insert/isaac_02.01/Grouper.1.4";
#my $vcfSuffix = "all.SVs";
my $vcfSuffix = "final.germline.SVs";
my @kidDirs  = qw(NA12881  NA12884  NA12887 
		  NA12879  NA12882  NA12885 
		  NA12888 
		  NA12880  NA12883  NA12886
		  );

my @parentDirs = qw(NA12877 NA12878 );

my @matGrandDirs = qw(NA12891 NA12892 );
my @patGrandDirs = qw(NA12889 NA12890 );

@matGrandDirs = ();
@patGrandDirs = ();

#@kidDirs = qw(NA12879 NA12880 NA12881 NA12882 );

my @allDirs = ( @kidDirs, @parentDirs, @matGrandDirs, @patGrandDirs );

#@allDirs = qw(NA12877 NA12878 NA12881);

## Set up 
foreach my $dir (@allDirs) {
    my $topDir  = "$buildDir/$dir";
    my $trioDir = "$buildDir/$dir/trioAnalysis.$rc";
    print STDERR "[Building]: $trioDir\n";
    mkdir($trioDir);

    my $cmd = "zcat $topDir/$vcfSuffix.vcf.gz | $binDir/vcfToBedSplit.pl -o $topDir/$vcfSuffix.converted";
    print STDERR "[Running]: $cmd\n";
    system($cmd);
}

foreach my $dir (@kidDirs) {
    my $kidDir  = "$buildDir/$dir";
    my $trioDir = "$buildDir/$dir/trioAnalysis.$rc";

    my $dadDir = "$buildDir/$parentDirs[0]";
    my $momDir = "$buildDir/$parentDirs[1]";

    my $cmd = "$binDir/windowBed -w 10 -a $kidDir/$vcfSuffix.converted.bed -b $dadDir/$vcfSuffix.converted.bed > $trioDir/intersect1.txt";
    print STDERR "[Running]: $cmd\n";
    system($cmd);
    $cmd = "$binDir/windowBed -w 10 -a $kidDir/$vcfSuffix.converted.bed -b $momDir/$vcfSuffix.converted.bed > $trioDir/intersect2.txt";
    print STDERR "[Running]: $cmd\n";
    system($cmd);

    $cmd = "$binDir/vcfCompare.pl -rc $rc -c $kidDir/$vcfSuffix.converted.vcf " .
	"-r $dadDir/$vcfSuffix.converted.vcf -i $trioDir/intersect1.txt " . 
	"-r $momDir/$vcfSuffix.converted.vcf -i $trioDir/intersect2.txt " .
	"-o $trioDir/results > $trioDir/results.txt ";
    print STDERR "[Running]: $cmd\n";
    system($cmd);
}

## End of file
