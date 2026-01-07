module BendersX

# Import dependencies
using JuMP
using MathOptInterface
using Printf
using LinearAlgebra
using SparseArrays
const MOI = MathOptInterface
using DataFrames
using Documenter

# Include source files
include("types.jl")
include("utils/utils.jl")
include("modules/modules.jl") 
include("problems/problems.jl") 

# Export abstract types
export AbstractBendersEnv
export AbstractData
export AbstractMaster
export AbstractOracle, AbstractOracleParam, AbstractTypicalOracle, BasicOracleParam
export AbstractBendersSeq, AbstractBendersBnB
export AbstractRootNodePreprocessing, root_node_processing!
export TerminationStatus, NotSolved, TimeLimit, Optimal, InfeasibleOrNumericalIssue
export TimeLimitException, UnexpectedModelStatusException, UndefError, AlgorithmException
export AbstractLoopState, AbstractLoopLog, AbstractLoopParam
export AbstractBendersSeqState, AbstractBendersSeqLog, AbstractBendersSeqParam
export AbstractBendersBnBState, AbstractBendersBnBLog, AbstractBendersBnBParam
export BendersSeqState, BendersSeqLog, BendersSeqParam
export BendersBnBState, BendersBnBLog, BendersBnBParam
export Hyperplane, aggregate, generate_cuts, set_parameter!, hyperplanes_to_expression, select_top_fraction, evaluate_violation, add_constraints
export get_sec_remaining, record_iteration!, update_upper_bound_and_gap!, is_terminated, check_lb_improvement!, print_iteration_info, to_dataframe

export AbstractNorm, StandardNorm, LpNorm
export SplitIndexSelectionRule, RandomFractional, MostFractional, LargestFractional
export DisjunctiveCutsAppendRule, NoDisjunctiveCuts, AllDisjunctiveCuts, DisjunctiveCutsSmallerIndices

end # module BendersX


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
