#! perl -w

my %count; 
my %sum;
open IN, "time.log";
open OUT, ">time.txt";
while (<IN>){
	if (/sim([A-Z])\/out(\d+)\/s(\d+)_([risq]+\d?)\.log\.txt.+=\s+([\d\.]+)\s+s/){
		my $group = $1 . $2 . $4;
		$count{$group} ++;
		$sum{$group} += $5;
		my $t = $5; 
		print OUT "$1\t$2\t$3\t$4\t$t\n"; 
	}	
}
close IN;
close OUT;

foreach my $g (qw/A B E F/){
	foreach my $k (1 .. 3){
		my $r = "$g$k" . 'r';
		my $a = "$g$k" . 'i1';
		my $b = "$g$k" . 'i2';
		my $c = "$g$k" . 'sq';
		printf "$g\t$k\t$count{$r}\t%.4f\t$count{$a}\t%.4f\t$count{$b}\t%.4f\t$count{$c}\t%.4f\n", $sum{$r}/$count{$r}, $sum{$a}/$count{$a}, $sum{$b}/$count{$b}, $sum{$c}/$count{$c};		
	}
}


