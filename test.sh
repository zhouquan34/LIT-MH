#! /bin/sh


./gimh -m test_data/test2.mat -p test_data/test.ph -w 0 -s 2000 -o test_out/try1 -r 34 --g-exp 3 --start 1 --kappa 1 -t 
./gimh -m test_data/test2.mat -p test_data/test.ph -s 500000 --omit test_data/test.skip -o test_out/try2 -r 34 --g-exp 3 --start 1 --kappa 1 -t  --long 0.05


