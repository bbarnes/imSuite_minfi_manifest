#! /usr/bin/perl -w

use strict;
use Getopt::Long;

use strict;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use File::Path qw(make_path remove_tree);
use File::Basename;
use PerlIO::gzip;
use Compress::Zlib;
use Storable;
use IPC::Open2;
use FileHandle;

use constant FALSE => 0;
use constant TRUE  => 1;

my $max_size = 48344810;
my $num_bins = 48;
my $out_dir  = "/Users/bretbarnes/Documents/data/improbe/scratch/cgnDB/dbSNP_Core4/design-output/improbe-cgn-top";

GetOptions("max=i" => \$max_size,
	   "bin=i" => \$num_bins,

	   "out=s" => \$out_dir,
    );

if (! defined($max_size) || ! defined($num_bins) || ! defined($out_dir)) {
     print STDERR "[Usage]: $0 -max max_size -bin num_bins -out out_dir [options] < STDIN\n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n\n";
    exit(1);
}

my $bin_size = int($max_size/$num_bins);

print STDERR "Will build $num_bins number of bins of size $bin_size!\n\n";

mkdir($out_dir) if (! -e $out_dir);
$out_dir = "$out_dir/bins-" . ($num_bins + 2);

mkdir($out_dir) if (! -e $out_dir);
print STDERR "Will write to out_dir: $out_dir\n\n";

# Clean all existing files...
system("rm $out_dir/*");

my @end = ();
my @beg = ();

my @fns = ();
my @fis = ();
my @fhs = ();

#
# TBD:: Write out report.tsv file
#
my $map_csv = "$out_dir/map.csv";
open(MAP, ">$map_csv") or die("Failed to open $map_csv for writing: $!");
print MAP "Index,Range_Beg,Range_End,Path\n";

for (my $ii = 0; $ii < $num_bins + 2; $ii++) {

    my $beg_pos = ($ii+0) * $bin_size + 1;
    my $end_pos = ($ii+1) * $bin_size;

    $beg[$ii] = $beg_pos;
    $end[$ii] = $end_pos;
    
    $fns[$ii] = "$out_dir/cgn2top.bin-$ii.range-$beg[$ii]-$end[$ii].tsv";
    $fis[$ii] = "$out_dir/cgn2top.bin-$ii.range-$beg[$ii]-$end[$ii].cgn-sorted.tsv.gz";

    # Clean all existing files...
    system("rm -f $fns[$ii]") if ( -e $fns[$ii] );
    system("rm -f $fis[$ii]") if ( -e $fis[$ii] );

    open($fhs[$ii], ">$fns[$ii]") or die("Failed to open $fns[$ii] for writing: $!");

    print MAP join(",", $ii,$beg[$ii],$end[$ii],$fis[$ii])."\n";
    
    print STDERR "[bin]:\tii=$ii; beg=$beg[$ii]; end=$end[$ii]; fns='$fns[$ii]'\n";
}
close(MAP);

my $header;
my $max_ten = int($max_size/100);
my $cur_cnt = 0;
my $max_cnt = 0;

while(<STDIN>) {
    s/^\s+//;
    s/\s+$//;
    next if /^(\s)*$/;
    # next if /^cgn/;

    my $line = $_;
    my @data = split(/\t/, $line);

    if ($line =~ m/^cgn/) {
	$header = $line;
	print STDERR "Setting header = '$header'...\n\n";
	next;
    }

    # cgn     epic    h450    h27k    h285    top
    my $cgn = $data[0];
    my $top = $data[5];

    my $idx = -1;
    for (my $ii = 0; $ii < @end; $ii++) {
	my $beg_pos = $beg[$ii];
	my $end_pos = $end[$ii];
	
	$idx = $ii if ($cgn >= $beg_pos && $cgn <= $end_pos);

	last if ($idx != -1);
    }
    if ($idx == -1) {
	print STDERR "[ERROR]: idx = $idx! cgn=$cgn, top=$top\n";
	print STDERR "[ERROR]:  line='$line'\n";
	exit(1);
    }

    print {$fhs[$idx]} join("\t", $cgn,$top)."\n";

    my $perc = int(100 * $cur_cnt / $max_size);
    print STDERR "[Status]: $perc ($cur_cnt); idx=$idx; cgn=$cgn, top=$top\n" if ( $cur_cnt % $max_ten == 0 );

    $cur_cnt++;
    last if ($max_cnt != 0 && $cur_cnt > $max_cnt);
}
# close($tarFH);

foreach my $fh (@fhs) {
    close($fh);
}
for (my $ii = 0; $ii < @fns; $ii++) {

    my $cmd = "sort -n -k 1,1 $fns[$ii] | gzip -c -> $fis[$ii]";
    print STDERR "[Running]: cmd='$cmd'\n";
    system($cmd);

    system("rm $fns[$ii]");
}
system("gzip $map_csv");

print STDERR "[Done]: Final Count = $cur_cnt\n\n";

exit(1);

# print STDERR "[Targets]: Loaded " . scalar(keys %tar) . "\n\n";
# print STDERR "[Loading]: Max Count=$maxCnt\n" if (defined($maxCnt));
# print STDERR "[Loading]: dataFile=$datFile\n";


#
# Subroutines::
#

sub trimSuffix {
    my $str=shift;

    for (my $ii = 0; $ii < 5; $ii++) {
	$$str=~ s/\.gz$//;  $$str=~ s/\.fa$//;
        $$str=~ s/\.bed$//; $$str=~ s/\.txt$//; $$str=~ s/\.csv$//; $$str=~ s/\.tsv$//;
        $$str=~ s/\.bam$//; $$str=~ s/\.sam$//; $$str=~ s/\.vcf$//; $$str=~ s/\.fas$//;
        $$str=~ s/\.bsp$//; $$str=~ s/\.manifest$//; $$str=~ s/\.manifests$//;
        $$str=~ s/\.sorted$//; $$str =~ s/\.merged$//; $$str =~ s/\.intersect$//;
	$$str=~ s/\.[^\.]-sorted$//;
    }
}

sub openFile {
    my $func=(caller(0))[3];
    my $fh=shift;
    my $file=shift;

    die("[$func]: ERROR File does not exist '$file'\n")
        if (! defined($file) || ! -e $file);

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    my $isBin = ($fileName =~ /\.bin$/) ? 1 : 0;
    trimSuffix(\$fileName);

    $$fh = gzopen($file, "rb") or die("Failed to open binary gzip stream=$file for reading: $!")
        if ($isBin && $fileSuffix eq ".gz");

    open($$fh, "<:gzip", $file) or die("Failed to open gzip stream=$file for reading: $!")
        if (! $isBin && $fileSuffix eq ".gz");

    open($$fh, "<:raw", $file) or die("Failed to open binary file=$file for reading: $!")
        if ($isBin && $fileSuffix ne ".gz");

    open($$fh, "<$file") or die("Failed to open regular file=$file for reading: $!")
	if (! $isBin && $fileSuffix ne ".gz");

    return(($fileName, $filePath, $fileSuffix, $isBin));
}


## End of file
