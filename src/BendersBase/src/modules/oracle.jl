export AbstractTypicalOracle, generate_cuts, set_parameter!, BasicOracleParam
"""
Abstract type for typical oracles used in Benders decomposition.
"""
abstract type AbstractTypicalOracle <: AbstractOracle end

"""
Prototype for the `generate_cuts` function.

Must be implemented by any concrete subtype of `AbstractOracle`. Given a candidate solution `(x_value, t_value)`, this method should attempt to separate the point via
valid inequalities.

Arguments:
- `x_value`: Given `x` solution.
- `t_value`: Given `t` solution.
- `tol_normalize`: Normalization tolerance for cut generation (default: 1.0).
- `time_limit`: Maximum time allowed for the oracle call (default: 3600 seconds).

Returns (to be implemented by concrete oracles):
- `is_in_L::Bool`: Whether the point is in the feasible region L (true if feasible, false if cuts were generated).
- `hyperplanes::Vector{Hyperplane}`: List of valid inequalities to be added to the master.
- `sub_obj_vals::Vector{Float64}`: Subproblem objective values for updating the upper bound. 
  Can be `NaN` if no meaningful objective was computed.

Throws an error if not implemented for a specific oracle type.
"""
function generate_cuts(oracle::AbstractTypicalOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    throw(UndefError("update generate_cuts for $(typeof(oracle))"))
end


"""
Basic parameter structure for oracles. Users can define oracle-specific parameter structures as subtypes of AbstractOracleParam. 
If the oracle has no specific parameter fields, use BasicOracleParam.
"""
struct BasicOracleParam <: AbstractOracleParam
    rtol::Float64
    atol::Float64
    zero_tol::Float64

    function BasicOracleParam(; rtol = 1e-9, atol = 0.0, zero_tol = 1e-9)
        new(rtol, atol, zero_tol)
    end
end

# Common utility functions for managing oracle parameters

function set_parameter!(oracle::AbstractOracle, param::AbstractOracleParam)
  if :oracle_param ∉ fieldnames(typeof(oracle))
      throw(UndefError("$(typeof(oracle)) must have a field named `oracle_param`"))
  elseif typeof(oracle.oracle_param) != typeof(param)
      throw(ArgumentError("Type mismatch: expected parameter of type $(typeof(oracle.oracle_param)), got $(typeof(param))"))
  else
      oracle.oracle_param = param
  end
end

function set_parameter!(oracle::AbstractOracle, param::String, value::Any)
  sym_param = Symbol(param)
  if sym_param ∈ fieldnames(typeof(oracle.oracle_param))
      setfield!(oracle.oracle_param, sym_param, value)
  else
      throw(ArgumentError("Parameter `$(param)` not found in `$(typeof(oracle.oracle_param))` for oracle of type `$(typeof(oracle))`"))
  end
end

include("oracleTypicalClassical.jl")
include("oracleTypicalSeparable.jl")
