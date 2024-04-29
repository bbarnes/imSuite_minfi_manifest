#! /usr/bin/perl -w

my $chunkSize = 1000000;

while(<STDIN>) {
    s/^\s+//;
    s/\s+$//;

    my $line = $_;
    my @data = split(/\t/, $line);

    my $chr   = $data[0];
    my $start = 1;
    my $end   = $data[1];

    for (my $ii = 1; $ii < $end; $ii += $chunkSize) {
	my $subEnd = $ii + $chunkSize;
	if ($subEnd >= $end) {
	    $subEnd = $end;
	    $ii = $end + 1;
	}
	$name = "$chr:$ii-$subEnd";
	print "$chr\t$ii\t$subEnd\t$name\t1.0\t+\n";
    }
}

## End of file
