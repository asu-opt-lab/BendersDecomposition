module BendersDecomposition

using Reexport

@reexport using BendersBase
@reexport using BendersLibrary


end # module BendersDecomposition

# To-Do: 
# 1. Remove mip.jl and Mip struct --> done
# 2. remove Data by having dim_x, dim_t, c_x, c_t info in master module --> done
# 3. clean initializers of master, oracle, env modules --> done
# 4. Kaiwen: refactor all files in script folder
# 5. Inho: refactor all snip-related files
# 6. rename `problem` to `data`
# 7. rename `oracle_param` to `param` for all oracles as we no longer has `solver_param`
# 8. rename `AbstractBendersDecomposition` to `AbstractBendersEnv`
# 9. Remove `AbstractCallbackParam` and `EmptyCallbackParam <: AbstractCallbackParam` and add `AbstractUserCallbackParam`; Lazy callback does not need parameters.
# 10. Return `to_dataframe(log)` for all solve! functions.
# 11. clean up all error-handling
# 12. Kaiwen: refactor root node preprocessing
# 13. Inho: unified oracle
# 14. Kaiwen: better handle infeasible subproblem for `ClassicalOracle` by incorporating normalization when infeasible
# 15. Inho: pareto oracle
