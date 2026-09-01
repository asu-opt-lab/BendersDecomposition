using JuMP
using CPLEX
using Gurobi

const COMPOSITION_SOLVERS = (:cplex, :gurobi)
const COMPOSITION_SUB_METHODS = (:default, :dual_simplex)

"""
    parse_solver_name(name) -> Symbol

Normalize and validate a solver name used by the solver-composition scripts.
"""
function parse_solver_name(name::AbstractString)
    solver = Symbol(lowercase(strip(name)))
    solver in COMPOSITION_SOLVERS || throw(
        ArgumentError(
            "Unsupported solver '$name'. Expected one of " *
            join(string.(COMPOSITION_SOLVERS), ", "),
        ),
    )
    return solver
end

"""
    parse_sub_method(name) -> Symbol

Normalize and validate the LP method used by subproblem optimizers.
"""
function parse_sub_method(name::AbstractString)
    method = Symbol(lowercase(strip(name)))
    method in COMPOSITION_SUB_METHODS || throw(
        ArgumentError(
            "Unsupported subproblem method '$name'. Expected one of " *
            join(string.(COMPOSITION_SUB_METHODS), ", "),
        ),
    )
    return method
end

"""
    composition_optimizer(solver, role; threads, seed, sub_method, tolerance)

Build a CPLEX or Gurobi optimizer for a solver-composition experiment.

Only resource limits, random seeds, and common numerical tolerances are fixed
for master optimizers. The master solver's branch-and-bound strategy remains at
its native default because it is the treatment being compared. For subproblem
optimizers, `sub_method = :dual_simplex` selects dual simplex in both solvers;
`:default` leaves the LP method at the vendor default.
"""
function composition_optimizer(
    solver::Symbol,
    role::Symbol;
    threads::Int,
    seed::Int,
    sub_method::Symbol,
    tolerance::Float64,
)
    solver in COMPOSITION_SOLVERS || throw(ArgumentError("Unsupported solver: $solver"))
    role in (:master, :subproblem) || throw(ArgumentError("Unsupported solver role: $role"))
    sub_method in COMPOSITION_SUB_METHODS || throw(ArgumentError("Unsupported subproblem method: $sub_method"))
    threads >= 1 || throw(ArgumentError("threads must be at least 1"))
    1e-9 <= tolerance <= 1e-2 || throw(ArgumentError("tolerance must lie in [1e-9, 1e-2]"))

    if solver == :cplex
        attributes = (
            "CPXPARAM_Threads" => threads,
            "CPXPARAM_RandomSeed" => seed,
            "CPX_PARAM_EPRHS" => tolerance,
            "CPX_PARAM_EPOPT" => tolerance,
            MOI.Silent() => true,
        )
        role == :master && (attributes = (attributes..., "CPX_PARAM_EPINT" => tolerance))
        role == :subproblem && sub_method == :dual_simplex &&
            (attributes = (attributes..., "CPXPARAM_LPMethod" => 2))
        return optimizer_with_attributes(CPLEX.Optimizer, attributes...)
    end

    attributes = (
        "Threads" => threads,
        "Seed" => seed,
        "FeasibilityTol" => tolerance,
        "OptimalityTol" => tolerance,
        MOI.Silent() => true,
    )
    role == :master && (attributes = (attributes..., "IntFeasTol" => tolerance))
    role == :subproblem && sub_method == :dual_simplex &&
        (attributes = (attributes..., "Method" => 1))
    return optimizer_with_attributes(Gurobi.Optimizer, attributes...)
end

function composition_optimizer(
    solver::AbstractString,
    role::Symbol;
    threads::Int,
    seed::Int,
    sub_method::AbstractString,
    tolerance::Float64,
)
    return composition_optimizer(
        parse_solver_name(solver),
        role;
        threads = threads,
        seed = seed,
        sub_method = parse_sub_method(sub_method),
        tolerance = tolerance,
    )
end
