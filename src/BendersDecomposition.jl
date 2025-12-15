module BendersDecomposition

using Reexport

@reexport using BendersBase
@reexport using BendersLibrary


end # module BendersDecomposition

# To-Do: 
# 4. Kaiwen: refactor all files in script folder
# 5. Inho: refactor all snip-related files
# 6. rename `problem` to `data` --> done
# 7. rename `oracle_param` to `param` for all oracles as we no longer has `solver_param` --> done
# 8. rename `AbstractBendersDecomposition` to `AbstractBendersEnv` --> done
# 9. rename `AbstractBendersCallback` to `AbstractBendersBnB`. --> done
# 10. rename `DisjunctiveOracle` to `SplitOracle`. --> done
# 11. rename the folders named `algorithms` to `envs` and include both `envs` and `oracles` folders inside `modules` folder. --> done
# 12. Remove `AbstractCallbackParam` and `EmptyCallbackParam <: AbstractCallbackParam` and add `AbstractUserCallbackParam`; Lazy callback does not need parameters. --> done
# 13. Return `to_dataframe(log)` for all solve! functions. --> done
# 14: Move output of disjunctive-cut statistics out of main code --> done
# 15. Refactor root node preprocessing --> done
# 16. clean up all error-handling
# 17. Inho: unified oracle
# 18. Kaiwen: better handle infeasible subproblem for `ClassicalOracle` by incorporating normalization when infeasible
# 19. Inho: pareto oracle
