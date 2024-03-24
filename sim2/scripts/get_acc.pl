#! perl -w

my @g = qw/r i1 i2 sq/;
foreach my $k (1 .. 20){
	foreach my $s (1 .. 5){
		foreach my $g (@g){
			print "$k\t$s\t$g";
			my $f = "dat/out/s${k}_${s}_$g.log.txt";
			open IN, "$f";
			while (<IN>){
				if (/^Add/){
					my @col = split /\s+/;
					foreach my $c (@col){
						$c =~ s/\(//;
						$c =~ s/\)//;
						$c =~ s/;//;
						if ($c =~ /\d/){
							print "\t$c";
						}
					}
				}
			}
			close IN;
			print "\n";
		}
	}
}



