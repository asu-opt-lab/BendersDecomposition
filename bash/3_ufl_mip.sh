#!/bin/sh
#SBATCH -t 0-01:00:00

echo "#!/bin/bash" > willdelete/jobscript_hm_$1.sh
# echo "#SBATCH -q grp_gbyeon" >> willdelete/jobscript_hm_$1.sh
# echo "#SBATCH -q gbyeonqos" >> willdelete/jobscript_hm_$1.sh
echo "#SBATCH -N 1" >> willdelete/jobscript_hm_$1.sh
echo "#SBATCH -n 20" >> willdelete/jobscript_hm_$1.sh
echo "#SBATCH -t 0-03:00:00" >> willdelete/jobscript_hm_$1.sh
echo "#SBATCH -o willdelete/hm.sh$1.out%j" >> willdelete/jobscript_hm_$1.sh
echo "#SBATCH -e willdelete/hm.sh$1.err%j" >> willdelete/jobscript_hm_$1.sh

echo "module purge" >>  willdelete/jobscript_hm_$1.sh
echo "module load julia" >> willdelete/jobscript_hm_$1.sh
echo "module load cplex" >> willdelete/jobscript_hm_$1.sh
echo "module load gurobi" >> willdelete/jobscript_hm_$1.sh

echo "julia --project=. bash/3_ufl_mip.jl --instance $1" >> willdelete/jobscript_hm_$1.sh

sbatch willdelete/jobscript_hm_$1.sh