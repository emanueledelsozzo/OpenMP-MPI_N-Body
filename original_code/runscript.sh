#!/bin/sh

make o3

./nbody > results.txt
./nbody >> results.txt
./nbody >> results.txt

python average.py results.txt
