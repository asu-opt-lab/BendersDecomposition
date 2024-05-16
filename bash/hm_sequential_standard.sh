for j in 1 2 3 5 7 10
    for ((i=1; i<=10; i++))
    do
        for cut_strategy in "ORDINARY_CUTSTRATEGY" "SPLIT_CUTSTRATEGY"
        do
            for SplitCGLPNormType in "STANDARDNORM" "L1GAMMANORM" "LINFGAMMANORM"
            do
                for SplitSetSelectionPolicy in "MOST_FRAC_INDEX" "RANDOM_INDEX"
                do
                    for StrengthenCutStrategy in "SPLIT_PURE_CUT_STRATEGY" "SPLIT_STRENGTHEN_CUT_STRATEGY" 
                    do
                        for SplitBendersStrategy in "NO_SPLIT_BENDERS_STRATEGY" "ALL_SPLIT_BENDERS_STRATEGY" "TIGHT_SPLIT_BENDERS_STRATEGY" 
                        do
                            bash bash/cfl_sequential.sh "f$(j)00-c(j)00-r5.0-p$i" $cut_strategy $SplitCGLPNormType $SplitSetSelectionPolicy $StrengthenCutStrategy $SplitBendersStrategy
                        done
                    done
                done
            done
        done
    done
done