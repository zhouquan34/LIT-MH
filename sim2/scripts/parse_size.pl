#! perl -w

my %count; 
my %sum;
open IN, "size.log";
while (<IN>){
	if (/dat\/out\/.+_(\d)_([risq]+\d?)\.log\.txt.+Mean\s+=\s+([\d\.]+);/){
		my $group = "s" . $1 . $2;
		$count{$group} ++;
		$sum{$group} += $3;
	}	
}
close IN;

foreach my $k (keys %count){
	printf "$k\t$count{$k}\t%.4f\n", $sum{$k}/$count{$k};
}


