
# ============================================================================
# Algorithm Strategy Types
# ============================================================================
# Solution Procedures
struct GenericSequential <: AbstractSequential 
    time_limit::Union{Float64,Nothing}
    iteration_limit::Union{Int,Nothing}
    gap_tolerance::Union{Float64}
    verbose::Bool
end

struct SequentialWithInout <: AbstractSequential 
    base::GenericSequential
    α::Float64
    λ::Float64
end

struct GenericCallback <: AbstractCallback 
    time_limit::Union{Float64,Nothing}
    gap_tolerance::Union{Float64}
    lazy_callback::Function
    user_callback::Union{Function, Nothing}
    verbose::Bool
end

struct SolverCallback{Solver} <: AbstractSolutionProcedure
    time_limit::Union{Float64,Nothing}
    gap_tolerance::Union{Float64}
    solver::Solver
    callback_function::Function
end

struct CallbackWithRootNodePreprocessing{T <: AbstractCallback, S <: AbstractSequential} <: AbstractCallback 
    base::T
    root_node_sequential::S
end

# TODO: it seems that we don't need this
# struct StochasticSequential <: AbstractSolutionProcedure end
# struct StochasticCallback <: AbstractSolutionProcedure end

# Cut Strategies
struct ClassicalCut <: AbstractCutStrategy end
# TODO: for the disjunctive benders, maybe we can have the performance restart?

########################################################
# Benders Problems
########################################################

# TODO: do not use type parameters, because the properties of the model are not known
# TODO: generative AI: because of multi dispatch
mutable struct GenericMasterProblem <: AbstractMasterProblem 
    model::Model
    variables::Dict{Symbol, Any}
    objective_value::Float64
    integer_variable_values::Vector{Float64}
    continuous_variable_values::Vector{Float64}
end

mutable struct GenericSubProblem <: AbstractSubProblem 
    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    other_constraints::Vector{ConstraintRef}
end

struct GenericMILP <: AbstractMILP end

########################################################
# Benders Environment
########################################################

mutable struct GenericBendersEnv <: AbstractBendersEnv
    master_problem::GenericMasterProblem
    sub_problem::GenericSubProblem
end






