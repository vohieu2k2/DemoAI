#!/bin/bash

wget http://snap.stanford.edu/data/web-Google.txt.gz
gunzip web-Google.txt.gz
../preproc web-Google.txt google.bin

data='./google.bin'
bin='../maxcut'
nodes='875713'
alg[1]='-A'
alg[2]='-E'
alg[3]='-M'
alg[4]='-Q'
alg[5]='-T'
alg[6]='-L'
nAlgs=6
#alg[1]="$1"
#nAlgs=1

# kminfrac='0.01'
# kmaxfrac='0.2'
# kincfrac='0.01'

# kmin=`echo "$nodes * $kminfrac" | bc`
# kmax=`echo "$nodes * $kmaxfrac" | bc`
# kinc=`echo "$nodes * $kincfrac" | bc`

kmin=10
kmax=210
kinc=50

echo "kmin=$kmin"
echo "kmax=$kmax"
echo "kinc=$kinc"
sleep 2

for k in `seq $kmin $kinc $kmax`;
do
    for a in `seq 1 1 $nAlgs`;
    do
	N='20'
	extra="-e 0.1 -r"
	F="-f"
	if [ "$a" -eq 2 ]; then
	    N='1'
	fi
	if [ "$a" -eq 5 ]; then
	    N='1'
	fi
	output="./resultsGoogle${alg[$a]}.txt"
	cmd="$bin -G $data ${alg[$a]} -o $output -k $k -N $N $extra $F"
	echo $cmd
	$cmd
    done
done

python3 plotGoogle.py   ## plot objective value
python3 plotGoogle.py q ## plot queries
python3 plotGoogle.py r ## plot rounds
