#!/bin/sh
#./testspme --model ../data/mono.mod --xyz ../data/N1000_Phi0.1.xyz \
#--dim 256 --porder 6 --xi 0.8 --rmax 8.371887

#./testspme --model ../data/mono.mod --xyz ../data/N1000_Phi0.1.xyz \
#--dim 256 --porder 8 --xi 0.8 --rmax 8.371887

#./testspme --model ../data/poly.mod --xyz ../data/N412.xyz \
#--dim 256 --porder 8 --xi 0.04 --rmax 200.37188

#./testspme --npts 1000 --xyz ./pvol01_40000.csv --dim 256 --porder 8 --xi 0.05 --rmax 300.0
#./testspme --npts 40000 --xyz ./pvol01_40000.csv --ref ./pvol01_40000_ref.csv --dim 256 --porder 8 --xi 0.04 --rmax 200.0
./testspme --npts 20000 --xyz ./pvol01_40000.csv --ref ./pvol01_20000_ref.csv --dim 256 --porder 8 --xi 0.04 --rmax 200.0
#./testspme --npts 5000 --xyz ./pvol01_40000.csv --ref ./pvol01_5000_ref.csv --dim 256 --porder 8 --xi 0.04 --rmax 200.0
