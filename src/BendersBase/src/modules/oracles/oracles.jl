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
- `tol_normalize`: Factor used to normalize the tolerance for cut generation (default: 1.0). This parameter is used only in the DCGLP procedure, and users typically do not need to modify it.
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
- `rtol` and `atol` determine the threshold for identifying a violation. Specifically, a violation is detected when
`sub_obj_val` >= `t_value` * (1 + `rtol`) + `atol`
in the finite-optimal subproblem case.
- `zero_tol` specifies the threshold below which a value is considered zero.
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

function set_parameter!(oracle::AbstractTypicalOracle, param::AbstractOracleParam)
  if :param ∉ fieldnames(typeof(oracle))
      throw(UndefError("$(typeof(oracle)) must have a field named `param`"))
  elseif typeof(oracle.param) != typeof(param)
      throw(ArgumentError("Type mismatch: expected parameter of type $(typeof(oracle.param)), got $(typeof(param))"))
  else
      oracle.param = param
  end
end

function set_parameter!(oracle::AbstractTypicalOracle, param::String, value::Any)
  sym_param = Symbol(param)
  if sym_param ∈ fieldnames(typeof(oracle.param))
      setfield!(oracle.param, sym_param, value)
  else
      throw(ArgumentError("Parameter `$(param)` not found in `$(typeof(oracle.param))` for oracle of type `$(typeof(oracle))`"))
  end
end

"""
    _validate_lp_compatibility(model::Model)

Validate that the model is compatible with typical oracles.
Typical oracles require a continuous Linear Programming (LP) subproblem:
- no integer/binary/semi-* variable-in-set constraints,
- an explicit linear objective,
- linear structural constraints with set type in {GreaterThan, LessThan, EqualTo},
- variable-in-set constraints limited to continuous bounds/fixes
  (GreaterThan, LessThan, EqualTo, Interval).

Throws UnsupportedModelException when unsupported model components are found.

This validation is shared by ClassicalOracle, UnifiedOracle, and ParetoOracle.
"""
function _validate_lp_compatibility(model::Model)
    # 1. Check for Integer/Binary variables (Mixed-Integer terms)
    for (_, S) in list_of_constraint_types(model)
        if S <: Union{MOI.Integer, MOI.ZeroOne, MOI.Semicontinuous, MOI.Semiinteger}
            throw(UnsupportedModelException(
                "Unsupported constraint type: $S. " *
                "Typical oracles require a continuous Linear Programming (LP) subproblem. " *
                "Integer or binary variables are not allowed."
            ))
        end
    end

    # 2. Check Objective Function (must be explicitly defined and linear)
    obj_type = objective_function_type(model)
    if obj_type === Nothing
        throw(UnsupportedModelException(
            "No objective function is defined. " *
            "Typical oracles require an explicit linear objective. " *
            "Use @objective(model, Min/Max, ...), e.g., @objective(model, Min, 0.0)."
        ))
    elseif !(obj_type <: Union{VariableRef, AffExpr, Real})
        throw(UnsupportedModelException(
            "Unsupported objective function type: $obj_type. " *
            "Typical oracles only support Linear Programming (LP) with linear objectives."
        ))
    end

    # 3. Check constraints in one pass (structural + variable-in-set)
    for con in all_constraints(model, include_variable_in_set_constraints=true)
        con_obj = constraint_object(con)
        set = con_obj.set
        func = con_obj.func

        # Variable-in-set constraints: allow continuous bounds/fixes only
        if func isa VariableRef
            if !(set isa MOI.GreaterThan || set isa MOI.LessThan || set isa MOI.EqualTo || set isa MOI.Interval)
                throw(UnsupportedModelException(
                    "Unsupported variable-in-set constraint type: $(typeof(set)). " *
                    "Typical oracles only support continuous variable bounds/fixes with " *
                    "GreaterThan (>=), LessThan (<=), EqualTo (==), or Interval."
                ))
            end
        # Structural constraints: linear function with supported scalar sets
        elseif func isa AffExpr
            if !(set isa MOI.GreaterThan || set isa MOI.LessThan || set isa MOI.EqualTo)
                throw(UnsupportedModelException(
                    "Unsupported constraint set type: $(typeof(set)). " *
                    "Typical oracles only support GreaterThan (>=), LessThan (<=), and EqualTo (==) constraints " *
                    "for affine structural constraints. If you have interval or other set types, please reformulate them."
                ))
            end
        else
            throw(UnsupportedModelException(
                "Unsupported constraint function type: $(typeof(func)). " *
                "Typical oracles only support linear constraints (LP). " *
                "Quadratic or other non-linear constraints are not supported."
            ))
        end
    end
end

include("oracleTypicalClassical.jl")
include("oracleTypicalSeparable.jl")
include("oracleTypicalUnified.jl")
include("oracleTypicalPareto.jl")
