#!/bin/sh
#SBATCH -t 0-01:00:00


bash bash/2_ufl_sequential_split.sh "KoerkelGhosh-asym/250/a/ga250a-1" 
bash bash/2_ufl_sequential_split.sh "KoerkelGhosh-asym/250/a/ga250a-2" 
bash bash/2_ufl_sequential_split.sh "KoerkelGhosh-asym/250/a/ga250a-3" 
bash bash/2_ufl_sequential_split.sh "KoerkelGhosh-asym/250/a/ga250a-4" 
bash bash/2_ufl_sequential_split.sh "KoerkelGhosh-asym/250/a/ga250a-5" 

