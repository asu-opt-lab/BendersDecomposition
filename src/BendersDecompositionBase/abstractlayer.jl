# ============================================================================
# Abstract type hierarchy
# ============================================================================

abstract type AbstractSolutionProcedure end
abstract type AbstractSequential <: AbstractSolutionProcedure end
abstract type AbstractCallback <: AbstractSolutionProcedure end
abstract type AbstractCutStrategy end


abstract type AbstractMasterProblem end
abstract type AbstractSubProblem end
abstract type AbstractMILP end
abstract type AbstractBendersEnv end

function solve!(env::AbstractBendersEnv, solution_procedure::AbstractSolutionProcedure, cut_strategy::AbstractCutStrategy)
    throw(MethodError(solve!, (env, solution_procedure, cut_strategy)))
end

function generate_cuts(env::AbstractBendersEnv, cut_strategy::AbstractCutStrategy, args...)
    throw(MethodError(_generate_cuts, (env, cut_strategy, args...)))
end





