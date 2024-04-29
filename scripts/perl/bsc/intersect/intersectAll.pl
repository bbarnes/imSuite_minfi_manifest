#! /usr/bin/perl -w

use strict;
use Getopt::Long;

my $launchExe = "/illumina/scratch/methylation/Polaris/shells/launchClusterLight.sh";
#my $launchExe = "/illumina/scratch/methylation/Polaris/shells/launchCluster.sh";

#my $intBedExe = "/illumina/thirdparty/BEDTools/BEDTools-Version-2.16.2/bin/intersectBed";
#my $intBedExe = "/illumina/thirdparty/from-uk-isilon/BEDTools/BEDTools-Version-2.17.0/bin/intersectBed";
my $intBedExe = "/usr/local/bin/intersectbed";

my $outDir;

my @regBeds = ();
my @desBeds = ();

my $verbose;
my $reverse;

##
## Subroutines
##

sub mssg($)
{
    my ($str) = @_;
    print STDERR $str if ($verbose);
}

sub buildDir($)
{
    my ($dir) = @_;

    if (-e $dir) {
        my $clean_cmd = "rm -rf $dir/*";
        mssg("\# Cleaning '$dir'\n");
        #system($clean_cmd);
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

##
## Main
##

GetOptions("o=s"   => \$outDir,
	   "r=s@"  => \@regBeds,
	   "d=s@"  => \@desBeds,
	   "v|verbose" => \$verbose,
	   "rev"       => \$reverse,
    );

if (scalar(@regBeds) == 0 || scalar(@desBeds) == 0 || ! defined ($outDir)) {
    print STDERR "[Usage]: $0 -o outDir -r regionBeds(s) -d desBed(s) [options]\n";
    print STDERR "[Usage]: Missing arguments! Exitng...\n";
    #print STDERR "@regBeds\n@desBeds\n$outDir\n";
    exit(1);
}

$outDir =~ s/\/+$//;
buildDir($outDir);

my $intCnt = 0;
foreach my $regBed (@regBeds) {
    #print STDERR "REG: $regBed\n";
    my $regName;
    if ($regBed =~ m/^.*\/([^\/]+)\.bed$/ || $regBed =~ m/^.*\/([^\/]+)\.bed.gz$/) {
	$regName = $1;
    } else {
	die("[ERROR]: Failed to parse region name from bed: $regBed!");
    }
    foreach my $desBed (@desBeds) {
	my ($desName, $shellFile);

	print STDERR "\# Ready to process design BED=$desBed\n";
	if ($desBed =~ m/^.*\/([^\/]+)\.bed$/ || $desBed =~ m/^.*\/([^\/]+)\.bed\.gz$/) {
	    $desName = $1;
	} else {
	    die("[ERROR]: Failed to parse design name from bed: $desBed!");
	}

	$shellFile = "$outDir/$desName.$regName.intersect.sh" if (! $reverse);;
	$shellFile = "$outDir/$regName.$desName.intersect.sh" if (  $reverse);;

	if (-e "$outDir/$desName.$regName.intersect.txt.gz") {
	    print STDERR "[Skipping]: $outDir/$desName.$regName.intersect.txt.gz already exists!\n";
	    next;
	}

	open(SHELL, ">$shellFile")
	    or die("Failded to open $shellFile for writing: $!");

	my $cmd;
	$cmd = "$intBedExe -loj -a $desBed -b $regBed > $outDir/$desName.$regName.intersect.txt" if (! $reverse);
	$cmd = "$intBedExe -loj -b $desBed -a $regBed > $outDir/$regName.$desName.intersect.txt" if (  $reverse);
	print SHELL "$cmd\n";
	$cmd = "gzip $outDir/$desName.$regName.intersect.txt";
	print SHELL "$cmd\n";
	close(SHELL);
	system("chmod 777 $shellFile");

	print STDERR "[CMD]: '$cmd'\n";
	system("$shellFile");
	#system("$launchExe int$intCnt $shellFile");
	$intCnt++;
	#exit(1);
    }
}

## End of file
