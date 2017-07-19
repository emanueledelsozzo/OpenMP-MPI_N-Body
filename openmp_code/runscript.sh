#!/bin/sh

make o3

threads=16

if [ "$#" == "1" ]; then
	threads=$1
fi

echo "num threads = $threads"

echo "warm up"

OMP_NUM_THREADS=$threads ./nbody > /dev/null
OMP_NUM_THREADS=$threads ./nbody > /dev/null

echo 0
OMP_NUM_THREADS=$threads ./nbody > results.txt
echo 1
OMP_NUM_THREADS=$threads ./nbody >> results.txt
echo 2
OMP_NUM_THREADS=$threads ./nbody >> results.txt
echo 3
OMP_NUM_THREADS=$threads ./nbody >> results.txt
echo 4
OMP_NUM_THREADS=$threads ./nbody >> results.txt

python average.py results.txt
diff goldenResults.dmp Population0009.dmp
diff goldenStage.jpg stage009.jpg
diff goldenPopulation.sta Population.sta
