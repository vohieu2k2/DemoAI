#!/bin/bash

../../ba 100 5 ba.100.5.txt
../../preproc ba.100.5.txt ba.100.5.bin
data='./ba.100.5.bin'
bin='../../maxcut'

alg[1]='-E'
nAlgs=1

nsamps[1]=0
nsamps[2]=1
nsamps[3]=10
nsamps[4]=100
#nsamps[5]=1000
nnsamps=4

k=10
#kmax=200
for a in `seq 1 1 $nnsamps`;
do
    N='1'
    extra="-e 0.1 -v -r"
    output="./resultsBA${alg[1]}-${nsamps[$a]}.txt"
    cmd="$bin -G $data ${alg[1]} -o $output -k $k -N $N $extra $F -n ${nsamps[$a]}"
    echo $cmd
    #$cmd
done

alg[1]='-A'
alg[2]='-F'
alg[3]='-T'
alg[4]='-L'
alg[5]='-M'
alg[6]='-Q'
nAlgs=6


kmin=10
#kmax=200
kmax=10
for k in `seq $kmin 10 $kmax`;
do
    for a in `seq 1 1 $nAlgs`;
    do
	N='20'
	extra="-e 0.1 -v -r"
	output="./resultsBA${alg[$a]}.txt"
	F='-f'
	cmd="$bin -G $data ${alg[$a]} -o $output -k $k -N $N $extra $F"
	echo $cmd
	$cmd
    done
done

python3 plotEne.py
python3 plotEne.py v
python3 plotEne.py q
