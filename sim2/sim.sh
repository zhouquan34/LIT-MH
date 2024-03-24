#! /bin/sh

k=$1	
sigma=(0.1 0.2 0.3 0.4 0.5)
dir="dat"	

for i in {1..5}
	do
	s=${sigma[$i-1]}
	Rscript sim_data.R $k $s 1000 5000 100 
	Rscript stepwise_true.R dat/cor${k}_$i
	../gimh -m $dir/cor${k}_$i.mat -p $dir/cor${k}_$i.ph --true $dir/cor${k}_$i.beta --init $dir/cor${k}_$i.im2 -r 34 --kappa 1 --g-exp 1 -t \
	   -o $dir/out/s${k}_${i}_r --rw  -s 200000  
	../gimh -m $dir/cor${k}_$i.mat -p $dir/cor${k}_$i.ph --true $dir/cor${k}_$i.beta --init $dir/cor${k}_$i.im2 -r 34 --kappa 1 --g-exp 1 -t \
	   -o $dir/out/s${k}_${i}_i2  -s 2000 --add-max 2 --del-max 1 --add-min -2 --del-min -2 
	../gimh -m $dir/cor${k}_$i.mat -p $dir/cor${k}_$i.ph --true $dir/cor${k}_$i.beta --init $dir/cor${k}_$i.im2 -r 34 --kappa 1 --g-exp 1 -t \
	   -o $dir/out/s${k}_${i}_i1  -s 2000 --add-max 1 --del-max 0 --add-min -1 --del-min -1 
	../gimh --lb -m $dir/cor${k}_$i.mat -p $dir/cor${k}_$i.ph --true $dir/cor${k}_$i.beta --init $dir/cor${k}_$i.im2 -r 34 --kappa 1 --g-exp 1 -t \
	   -o $dir/out/s${k}_${i}_sq  -s 2000 
	rm dat/cor${k}_${i}.mat
done


