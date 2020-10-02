#!/bin/bash

make
./preproc ./data/ba.1k.10.txt ./data/ba.1k.10.bin

echo "generating table in README..."
echo 

data='./data/ba.1k.10.bin'
bin='./maxcut'
k='200'
epsi='0.1'
repetitions='100'

printf "%20s %20s %20s %20s\n" "Algorithm" "Objective" "Queries" "Rounds"

#run IteratedGreedy
out=`$bin -G $data -k $k -q -e $epsi -T -f -r`
printf "%20s %20.0f %20.0f %20.0f\n" "IteratedGreedy" $out

#run Ene et al. 2020 (with exact oracles for multilinear extension)
out=`$bin -G $data -k $k -q -e $epsi -f -r -E`
printf "%20s %20.0f %20.0f %20.0f\n" "Ene et al. 2020" $out

#run FastRandomGreedy
out=`$bin -G $data -k $k -q -e $epsi -f -r -Q -N $repetitions`
printf "%20s %20.0f %20.0f %20.0f\n" "FastRandomGreedy" $out

#run AdaptiveNonomonotoneMax
out=`$bin -G $data -k $k -q -e $epsi -f -r -M -N $repetitions`
printf "%20s %20.0f %20.0f %20.0f\n" "AdaNonmonotoneMax" $out

#run AdaptiveSimpleThreshold
out=`$bin -G $data -k $k -q -e $epsi -f -r -A -N $repetitions`
printf "%20s %20.0f %20.0f %20.0f\n" "AdaSimpleThresh" $out

#run AdaptiveThresholdGreedy
out=`$bin -G $data -k $k -q -e $epsi -f -r -L -N $repetitions`
printf "%20s %20.0f %20.0f %20.0f\n" "AdaThreshGreedy" $out
