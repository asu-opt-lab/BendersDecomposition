#!/bin/sh
#SBATCH -t 0-01:00:00


for ((i=1; i<=10; i++))
do
    bash bash/1_cfl_sequential_standard.sh "f700-c700-r5.0-p$i" 
done
