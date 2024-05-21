#!/bin/sh
#SBATCH -t 0-01:00:00


for ((i=1; i<=2; i++))
do
    for cut_strategy in "SPLIT_CUTSTRATEGY"
    do
        for SplitCGLPNormType in "STANDARDNORM" "L1GAMMANORM" "LINFGAMMANORM"
        do
            for SplitSetSelectionPolicy in "MOST_FRAC_INDEX" "RANDOM_INDEX"
            do
                for StrengthenCutStrategy in "SPLIT_PURE_CUT_STRATEGY" "SPLIT_STRENGTHEN_CUT_STRATEGY" 
                do
                    for SplitBendersStrategy in "NO_SPLIT_BENDERS_STRATEGY" "ALL_SPLIT_BENDERS_STRATEGY" "TIGHT_SPLIT_BENDERS_STRATEGY" 
                    do
                        bash bash/cfl_sequential_split.sh "f500-c500-r5.0-p$i" $cut_strategy $SplitCGLPNormType $SplitSetSelectionPolicy $StrengthenCutStrategy $SplitBendersStrategy
                    done
                done
            done
        done
    done
done
