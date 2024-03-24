#! /bin/bash


dir=$1
pre=$2
n=$3
p=$4
k=$5

snrs=(3 2 1)

for m in {1..3}
	do
	snr=${snrs[$m-1]}
	echo $m $snr
	Rscript sim_data_${pre}.R $k $snr $n $p $dir $m
	../gimh -m $dir/${pre}${k}_$m.mat -p $dir/${pre}${k}_${m}.ph --true $dir/${pre}${k}_${m}.beta -r $k  --kappa 2 --g-exp 3 -t --start 10 -s 100000 -o $dir/out$m/s${k}_r  --rw 
	../gimh -m $dir/${pre}${k}_$m.mat -p $dir/${pre}${k}_${m}.ph --true $dir/${pre}${k}_${m}.beta -r $k  --kappa 2 --g-exp 3 -t --start 10 -s 2000   -o $dir/out$m/s${k}_i1  --add-max 1 --del-max 0 --add-min -1 --del-min -1  
	../gimh -m $dir/${pre}${k}_$m.mat -p $dir/${pre}${k}_${m}.ph --true $dir/${pre}${k}_${m}.beta -r $k  --kappa 2 --g-exp 3 -t --start 10 -s 2000   -o $dir/out$m/s${k}_i2  --add-max 2 --del-max 1 --add-min -2 --del-min -2  
	../gimh --lb -m $dir/${pre}${k}_$m.mat -p $dir/${pre}${k}_${m}.ph --true $dir/${pre}${k}_${m}.beta -r $k  --kappa 2 --g-exp 3 -t --start 10 -s 2000   -o $dir/out$m/s${k}_sq   
	rm $dir/${pre}${k}_${m}.*
done

