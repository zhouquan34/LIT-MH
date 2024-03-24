#! /bin/sh


mkdir jobs

for g in A B E F
	do
	mkdir sim$g
	for j in 1 2 3
		do
		mkdir sim$g/out$j
	done
done

