#!/bin/sh

#make o3

threads=16
proc=1

if [ "$#" == "1" ]; then
        proc=$1
fi

if [ "$#" == "2" ]; then
	proc=$1
	threads=$2
fi

echo "num processes = $proc"
echo "num threads = $threads"

echo "warm up"

mpirun -np $proc -x OMP_NUM_THREADS=$threads ./nbody > /dev/null
mpirun -np $proc -x OMP_NUM_THREADS=$threads ./nbody > /dev/null

echo 0
mpirun -np $proc -x OMP_NUM_THREADS=$threads ./nbody > results.txt
echo 1
mpirun -np $proc -x OMP_NUM_THREADS=$threads ./nbody >> results.txt
echo 2
mpirun -np $proc -x OMP_NUM_THREADS=$threads ./nbody >> results.txt
echo 3
mpirun -np $proc -x OMP_NUM_THREADS=$threads ./nbody >> results.txt
echo 4
mpirun -np $proc -x OMP_NUM_THREADS=$threads ./nbody >> results.txt

echo "num processes = $proc" >> final_results.txt
echo "num threads = $threads" >> final_results.txt

python average.py results.txt >> final_results.txt
diff goldenResults.dmp Population0009.dmp
diff goldenStage.jpg stage009.jpg
diff goldenPopulation.sta Population.sta
