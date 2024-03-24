#! perl -w

my %map = (
	'r' => 0,
	'i1' => 1,
	'i2' => 2,
	'sq' => 3
);

our @stats;

sub read_table{
	my $f = $_[0];
	open IN, $f;
	while (<IN>){
		chomp;
		my @col = split /\s+/;
		if (exists $map{$col[0]}){
			my $g = $map{$col[0]};
			foreach my $k (2 .. $#col){
				push @{$stats[$g]->[$col[1] - 1]}, $col[$k]; 
			}
		}
	}
	close IN; 
}

read_table("time.table");

open IN, "scripts/modes.txt";
while (<IN>){
	chomp;
	my @col = split /\s+/;
	if (exists $map{$col[0]}){
		my $g = $map{$col[0]};
		foreach my $i (1..5){
			push @{$stats[$g]->[$i-1]}, $col[$i];
		}		
	}
}
close IN;

read_table("acc.table");
read_table("ess.table");
read_table("size.table");

open OUT, ">stats2.tex";

print OUT <<HERE; 

\\begin{tabular}{cccccc}
\\toprule
  &  &  Random walk & LIT1 & LIT2 & LIB  \\\\
\\midrule
HERE


my @names = ('Time', 'Local modes', 'Acc. Rate',  'ESS($T_1$)/Time',  'ESS($T_2$)/Time');
my @sigma = qw/0.1 0.2 0.3 0.4 0.5/;
foreach my $nr (0 .. 24){
	my $g2 = int($nr / 5); 
	my $r2 = $nr % 5;  
	if ($r2 == 0){
		if ($g2 > 0){print OUT "\\cmidrule{1-6} \n"; }
		print OUT ' \\multirow{5}{*}{ \\shortstack{ $\\sigma_\beta = ' . $sigma[$g2] . 
		  '$ \\\\ Mean model size = ', sprintf("%.1f", $stats[1]->[$g2]->[5]), ' }';
	}else{
		print OUT ' ';
	}
	print OUT " & $names[$r2] ";	
	foreach my $k (0 .. 3){
		if ($r2 <= 2){
			print OUT ' & ', formatn($stats[$k]->[$g2]->[$r2]);
		}else{
			print OUT ' & ', formatn($stats[$k]->[$g2]->[$r2] / $stats[$k]->[$g2]->[0]);
		}
	}
	print OUT "\\\\ \n";
}	

print OUT <<HERE; 
\\bottomrule
\\end{tabular}
HERE


close OUT;

sub formatn{
	my $x = $_[0];
	if ($x > 10){
		return sprintf("%.1f", $x);
	}elsif ($x > 1){
		return sprintf("%.2f", $x);
	}elsif ($x > 0.01){
		return sprintf("%.3f", $x);
	}else{
		return sprintf("%.4f", $x);
	}
}


