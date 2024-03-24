#! perl -w

my $nk=3;
my @snr = qw/3 2 1/;
my @data = (
	"\$n = 500, p = 1000\$, \\\\ independent design", 
	"\$n = 500, p = 1000\$, \\\\ correlated design", 
	"\$n = 1000, p = 5000\$, \\\\ independent design", 
	"\$n = 1000, p = 5000\$, \\\\ correlated design"
	);
my %map = ('A' => 0,
	'B' => 1,
	'E' => 2,
	'F' => 3);

my @stats; 
my @mods;

open IN, "time.table";
while (<IN>){
	chomp;
	my @col = split /\s+/;
	my $g = $map{$col[0]};
	foreach my $k (0 .. $nk){
		push @{$stats[$g]->[$col[1] - 1]->[$k]}, $col[2*$k + 3]; 
	}
}
close IN; 

open IN, "fail.table";
while (<IN>){
	chomp;
	my @col = split /\s+/;
	my $g = $map{$col[0]};
	foreach my $k (0 .. $nk){
		push @{$stats[$g]->[$col[1] - 1]->[$k]}, 100 - $col[$k + 2]; 
	}
}
close IN; 


open IN, "hit.table";
while (<IN>){
	chomp;
	my @col = split /\t/;
	my $g = $map{$col[0]};
	foreach my $k (0 .. $nk){
		push @{$stats[$g]->[$col[1] - 1]->[$k]}, $col[$k + 4]; 
	}
	my @c = @col[2 .. 3];
	$mods[$g]->[$col[1]-1] = \@c;
}
close IN; 

open IN, "Tmax.txt";
while (<IN>){
	chomp;
	my @col = split /\s+/;
	my $g = $map{$col[0]};
	foreach my $k (0 .. $nk){
		push @{$stats[$g]->[$col[1] - 1]->[$k]}, $col[$k + 2]; 
	}	
}
close IN;



open OUT, ">stats.tex";

print OUT <<HERE; 

\\begin{tabular}{p{3.5cm}cccccc}
\\toprule
 &  &  &  Random walk & LIT1 & LIT2 & Bal \\\\
\\midrule
HERE


my @names = ("Time",  "Success",  '$H_{\\rm{max}}$', '$t_{\\rm{max}}$');
foreach my $nr (0 .. 47){
	my $g1 = int($nr / 12); # ABEF
	my $r1 = $nr % 12;  
	my $g2 = int($r1 / 4); # snr,123
	my $r2 = $r1 % 4;  
	if ($r1 == 0){
		if ($g1 > 0){print OUT "\\midrule \n"; }
		print OUT '\\multirow{12}{*}{\\shortstack{' . $data[$g1] . '}}';
	}else{
		print OUT ' '; 
	}	
	if ($r2 == 0){
		if ($g2 > 0){print OUT "\\cmidrule{2-7} \n"; }
		print OUT ' & \\multirow{4}{*}{ SNR = ' . $snr[$g2] . '} ';
	}else{
		print OUT ' & ';
	}
	print OUT " & $names[$r2] ";	
	foreach my $k (0 .. $nk){
		print OUT ' & ', $stats[$g1]->[$g2]->[$k]->[$r2];
	}
	print OUT "\\\\ \n";
}	

print OUT <<HERE; 
\\bottomrule
\\end{tabular}
HERE


close OUT;


