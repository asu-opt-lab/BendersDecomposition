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
# 8. rename `AbstractBendersDecomposition` to `AbstractBendersEnv`
# 9. rename \texttt{AbstractBendersCallback} to \texttt{AbstractBendersBnB}.
# 10. rename \texttt{DisjunctiveOracle} to \texttt{SplitOracle}.
# 11. rename the folder \texttt{algorithms} to \texttt{envs}.
# 12. Remove `AbstractCallbackParam` and `EmptyCallbackParam <: AbstractCallbackParam` and add `AbstractUserCallbackParam`; Lazy callback does not need parameters.
# 13. Return `to_dataframe(log)` for all solve! functions.
# 14. clean up all error-handling
# 15. Kaiwen: refactor root node preprocessing
# 16. Inho: unified oracle
# 17. Kaiwen: better handle infeasible subproblem for `ClassicalOracle` by incorporating normalization when infeasible
# 18. Inho: pareto oracle
