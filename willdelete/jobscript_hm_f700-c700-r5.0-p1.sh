#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 0-03:00:00
#SBATCH -o willdelete/hm.shf700-c700-r5.0-p1.out%j
#SBATCH -e willdelete/hm.shf700-c700-r5.0-p1.err%j
module purge
module load julia
module load cplex
module load gurobi
julia --project=. bash/0_cfl_sequential_split.jl --instance f700-c700-r5.0-p1
