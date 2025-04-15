export AbstractTypicalOracle, AbstractDisjunctiveOracle, generate_cuts

"""
Abstract type for typical oracles used in Benders decomposition.
"""
abstract type AbstractTypicalOracle <: AbstractOracle end

"""
Abstract type for disjunctive oracles
"""
abstract type AbstractDisjunctiveOracle <: AbstractOracle end

"""
Prototype for the `generate_cuts` function.

Must be implemented by any concrete subtype of `AbstractOracle`. Given a candidate solution `(x_value, t_value)`, this method should attempt to separate the point via
valid inequalities.

Arguments:
- `x_value`: Given `x` solution.
- `t_value`: Given `t` solution.
- `tol`: Numerical tolerance for cut generation (default: `1e-6`).
- `time_limit`: Maximum time allowed for the oracle call (default: 3600 seconds).

Returns (to be implemented by concrete oracles):
- `is_separated::Bool`: Whether any cuts were generated.
- `hyperplanes::Vector{Hyperplane}`: List of valid inequalities to be added to the master.
- `sub_obj_vals::Vector{Float64}`: Subproblem objective values for updating the upper bound. 
  Can be `NaN` if no meaningful objective was computed.

Throws an error if not implemented for a specific oracle type.
"""
function generate_cuts(oracle::AbstractOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol = 1e-6, time_limit = 3600)
    throw(UndefError("update generate_cuts for $(typeof(AbstractOracle))"))
end

include("oracleTypicalClassical.jl")
include("oracleTypicalSeparable.jl")
include("oracleDisjunctive.jl")