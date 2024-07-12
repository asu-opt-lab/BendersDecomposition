#!/bin/sh
#SBATCH -t 0-01:00:00


for ((i=1; i<=10; i++))
do
    for cut_strategy in "SPLIT_CUTSTRATEGY"
    do
        bash bash/0_cfl_sequential_split.sh "f200-c200-r5.0-p$i" 
    done
done
