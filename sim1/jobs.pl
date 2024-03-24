#! perl -w

my @dir = qw/simA simB simE simF/;
my @pre = qw/ind cor ind cor/;
my @ns = qw/500 500 1000 1000/;
my @ps = qw/1000 1000 5000 5000/;

foreach my $i (0 .. 3){
	my $dir = $dir[$i];
	my $pre = $pre[$i];
	my $n = $ns[$i];
	my $p = $ps[$i]; 
	foreach my $j (1 .. 5){
		open OUT, ">jobs/s${i}_$j.sh";
		print OUT  <<HERE; 
#! /bin/bash

HERE
		
		foreach my $m (1 .. 20){
			my $k = ($j-1)*20 + $m; 
			my $c = $i + 1;
			print OUT "bash simulation.sh $dir $pre $n $p $k \n"; 
		}

		close OUT;
		
		system("nohup bash jobs/s${i}_$j.sh &");	
	}
}



