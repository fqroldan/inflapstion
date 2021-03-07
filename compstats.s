#!/bin/bash

N=$1
param=$2

sender(){
	mkdir -p ../Output/CompStats/$param/run$1
	cp comp_stats.jl ../Output/CompStats/$param/run$1/comp_stats.jl
	cd ../Output/CompStats/$param/run$1
	mkdir -p ../Graphs/tests/
	sleep 1
	echo "started run $1"
	export JULIA_NUM_THREADS=20
	julia comp_stats.jl $1 $2 $3
	cd ../../../../Codes
}

rm -rf ../Output/CompStats/$param

for ((i=1; i<=$N; i++))
do
	sender $i $N $param
done

wait

echo "All done"