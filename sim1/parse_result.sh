#! /bin/bash

#for s in 'A' 'B' 'E' 'F'
#	do
#	for i in {1..3}
#		do
#		echo $s $i
#		Rscript scripts/first_hit.R sim$s/out$i/s >sim$s/out${i}.txt
#		Rscript scripts/average.R $s $i
#	done
#done

#grep "Time used" sim*/out*/*log* >time.log
#perl scripts/parse_time.pl | sort -k1n >time.table

#Rscript scripts/hit_table.R >hit.table
#Rscript scripts/compare_hit.R >fail.table
#Rscript scripts/make_table_t.R >Tmax.txt
#rm time.log time.txt

