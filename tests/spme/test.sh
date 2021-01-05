#!/bin/sh
#./testspme --model ../data/mono.mod --xyz ../data/N1000_Phi0.1.xyz \
#--dim 256 --porder 6 --xi 0.8 --rmax 8.371887

#./testspme --model ../data/mono.mod --xyz ../data/N1000_Phi0.1.xyz \
#--dim 256 --porder 8 --xi 0.8 --rmax 8.371887

#./testspme --model ../data/poly.mod --xyz ../data/N412.xyz \
#--dim 256 --porder 8 --xi 0.04 --rmax 200.37188

#./testspme --npts 1000 --xyz ./pvol01_40000.csv --dim 256 --porder 8 --xi 0.05 --rmax 300.0
#./testspme --npts 40000 --xyz ./pvol01_40000.csv --ref ./pvol01_40000_ref_scaled.csv --dim 256 --porder 8 --xi 0.04 --rmax 200.0
#./testspme --npts 20000 --xyz ./pvol01_40000.csv --ref ./pvol01_20000_ref_scaled.csv --dim 256 --porder 8 --xi 0.04 --rmax 200.0
#./testspme --npts 5000 --xyz ./pvol01_40000.csv --ref ./pvol01_5000_ref_scaled.csv --dim 256 --porder 8 --xi 0.04 --rmax 200.0



#./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 64  --porder 8 --xi 0.05 --rmax 100.0
#./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 128 --porder 8 --xi 0.05 --rmax 100.0
#./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.05 --rmax 100.0
#./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 8 --xi 0.05 --rmax 100.0


./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 6  --xi 0.05 --rmax 100.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 8  --xi 0.05 --rmax 100.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 10 --xi 0.05 --rmax 100.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 12 --xi 0.05 --rmax 100.0

./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 150.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 140.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 130.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 120.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 110.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 100.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 90.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 80.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 70.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 60.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 50.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 40.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 30.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 20.0
./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 256 --porder 8 --xi 0.03 --rmax 10.0





#./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 8 --xi 0.09 --rmax 200.0
#./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 8 --xi 0.09 --rmax 180.0
#./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 8 --xi 0.09 --rmax 160.0
#./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 8 --xi 0.09 --rmax 140.0
#./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 8 --xi 0.09 --rmax 120.0
#./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 8 --xi 0.09 --rmax 100.0
#./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 8 --xi 0.09 --rmax 80.0
#./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 8 --xi 0.09 --rmax 60.0
#./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 8 --xi 0.09 --rmax 40.0
#./testspme2 --npts 160000 --xyz ./data/vol01/pvol01_160000.csv --ref ./data/vol01_ref/pvol01_160000_ref_scaled.csv --dim 512 --porder 8 --xi 0.09 --rmax 20.0
