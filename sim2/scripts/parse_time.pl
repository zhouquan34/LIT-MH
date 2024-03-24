#! perl -w

my %count; 
my %sum;
open IN, "time.log";
while (<IN>){
	if (/dat\/out\/.+_(\d)_([risq]+\d?)\.log\.txt.+=\s+([\d\.]+)\s+s/){
		my $group = "s" . $1 . $2;
		$count{$group} ++;
		$sum{$group} += $3;
	}	
}
close IN;

foreach my $k (keys %count){
	printf "$k\t$count{$k}\t%.4f\n", $sum{$k}/$count{$k};
}


