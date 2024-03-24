#! perl -w

my %count;
my %acc; 
open IN, "acc.txt";
while (<IN>){
	my @col = split /\s+/;
	my $s = $col[1];
	my $g = $col[2];
	my @n = @col[(4, 6, 8)];
	my @acc = @col[3, 5, 7];
	my $sum = sum(@n);
	my $total = 2000;
	if ($g eq 'r'){$total = 200000;}	
	if ($total != $sum){
		print "$s\t$g\t$sum\n";
	}
	my $na = 0;
	map {$na += $n[$_] * $acc[$_] }0 .. 2;
	$count{$g}->[$s - 1] += $sum;
	$acc{$g}->[$s - 1] += $na;
}
close IN;

sub sum{
	my @x = @_;
	my $s = 0;
	map {$s += $_} @x;
	return $s
}

open OUT, ">acc.table";
foreach my $k (keys %count){
	foreach my $s (0 .. 4){
		my $ss = $s + 1;
		printf OUT "$k\t$ss\t%.4f\n", $acc{$k}->[$s] / $count{$k}->[$s];
	}
}
close OUT; 

open IN, "time_ave.txt";
open OUT, ">time.table";
while (<IN>){
	my @col = split /\s+/;
	if ($col[0] =~ /s(\d+)(.+)$/ ){
		my $s = $1;
		my $g = $2;
		print OUT "$g\t$s\t$col[2]\n";
	}	
}
close IN;
close OUT;

open IN, "size.txt";
open OUT, ">size.table";
while (<IN>){
	my @col = split /\s+/;
	if ($col[0] =~ /s(\d+)(.+)$/ ){
		my $s = $1;
		my $g = $2;
		print OUT "$g\t$s\t$col[2]\n";
	}	
}
close IN;
close OUT;


