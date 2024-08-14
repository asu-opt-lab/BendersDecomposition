#!/bin/sh
#SBATCH -t 0-01:00:00


for ((i=1; i<=10; i++))
do
    bash bash/0_cfl_sequential_split.sh "f1000-c1000-r5.0-p$i" 
done
