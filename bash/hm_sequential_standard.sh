#!/bin/sh
#SBATCH -t 0-01:00:00


for ((i=1; i<=2; i++))
do
    bash bash/cfl_sequential_standard.sh "f500-c500-r5.0-p$i" 
done
