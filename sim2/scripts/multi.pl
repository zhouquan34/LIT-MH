#! perl -w

my $g = $ARGV[0];
my @count = qw/0 0 0 0 0/;

foreach my $i (1 .. 20){
#foreach my $i (1 .. 1){
	foreach my $k (1 .. 5){
		open IN, "dat/out/s${i}_${k}_$g.path.txt";
		my %modes;
		while (<IN>){
			if ($. == 1){next;}
			chomp;
			my @col = split /\s+/;
			if ($col[$#col] == 1){
				$modes{$col[4]} = 1
			}
		}
		close IN;
		my $n = scalar(keys %modes);
		$count[$k-1] += $n;
		#print "$i\t$k\t$n\t", join(",", keys %modes), "\n";
	}
}

map {$count[$_] /= 20.0 } 0 .. 4;

print join("\t", @count), "\n";


