export AbstractTypicalOracle, generate_cuts, set_parameter!, BasicOracleParam, _validate_lp_compatibility, _normalize_to_scalar_constraints!
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
- no discontinuous(e.g., integer/binary/semi-*) variables,
- affine objective and constraints,

Throws `UnsupportedModelException` when unsupported model components are found.

This validation is shared by `ClassicalOracle`, `UnifiedOracle`, and `ParetoOracle`.
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

    # 2. Check Objective Function (must be affine)
    obj_type = objective_function_type(model)
    if !(obj_type <: Union{VariableRef, AffExpr, Real})
        throw(UnsupportedModelException(
            "Unsupported objective function type: $obj_type. " *
            "Typical oracles only support Linear Programming (LP) with linear objectives."
        ))
    end

    # 3. Check constraints
    for con in all_constraints(model, include_variable_in_set_constraints=true)
        con_obj = constraint_object(con)
        set = con_obj.set
        func = con_obj.func

        if func isa VariableRef || func isa AffExpr || (func isa AbstractVector && all(x -> x isa AffExpr || x isa VariableRef, func))
            if !(set isa MOI.GreaterThan || set isa MOI.LessThan || set isa MOI.EqualTo || set isa MOI.Interval || set isa MOI.Zeros || set isa MOI.Nonpositives || set isa MOI.Nonnegatives)
                throw(UnsupportedModelException(
                    "Unsupported constraint set type: $(typeof(set)). " *
                    "Typical oracles only support affine constraints."
                ))
            end
        else
            throw(UnsupportedModelException(
                "Unsupported constraint function type: $(typeof(func)). " *
                "Typical oracles only support affine constraints."
            ))
        end
    end
end

"""
    _insert_suffix(name::String, suffix::String) -> String

Insert suffix before the index brackets in a constraint name.
Examples:
- `demand[1]` + `_lb` -> `demand_lb[1]`
- `flow[i,j]` + `_ub` -> `flow_ub[i,j]`
- `simple` + `_lb` -> `simple_lb`
"""
function _insert_suffix(name::String, suffix::String)
    bracket_pos = findfirst('[', name)
    if bracket_pos === nothing
        return name * suffix
    else
        return name[1:bracket_pos-1] * suffix * name[bracket_pos:end]
    end
end

"""
    _normalize_to_scalar_constraints!(model::Model)

Pre-process all non-scalar constraints into scalar form so that oracle
transformations only encounter `GreaterThan`, `LessThan`, and `EqualTo`.

Splits:
- `MOI.Interval` → `GreaterThan(lower)` + `LessThan(upper)`
- `MOI.Zeros` → individual `EqualTo(0)` per row
- `MOI.Nonnegatives` → individual `GreaterThan(0)` per row
- `MOI.Nonpositives` → individual `LessThan(0)` per row
"""
function _normalize_to_scalar_constraints!(model::Model)
    # Collect constraints to split (avoid modifying during iteration)
    cons_to_split = ConstraintRef[]
    for con in all_constraints(model, include_variable_in_set_constraints=false)
        set = constraint_object(con).set
        if set isa MOI.Interval || set isa MOI.Zeros || set isa MOI.Nonnegatives || set isa MOI.Nonpositives
            push!(cons_to_split, con)
        end
    end

    for con in cons_to_split
        con_obj = constraint_object(con)
        set = con_obj.set
        func = con_obj.func
        original_name = name(con)

        if set isa MOI.Interval
            # lb <= expr <= ub  →  expr >= lb  AND  expr <= ub
            lb_con = @constraint(model, func >= set.lower)
            ub_con = @constraint(model, func <= set.upper)
            if !isempty(original_name)
                set_name(lb_con, _insert_suffix(original_name, "_lb"))
                set_name(ub_con, _insert_suffix(original_name, "_ub"))
            end

        elseif set isa MOI.Zeros
            # [f1; f2; ...] ∈ Zeros  →  f1 == 0, f2 == 0, ...
            for (i, fi) in enumerate(func)
                new_con = @constraint(model, fi == 0)
                if !isempty(original_name)
                    set_name(new_con, _insert_suffix(original_name, "_$i"))
                end
            end

        elseif set isa MOI.Nonnegatives
            # [f1; f2; ...] ∈ Nonneg  →  f1 >= 0, f2 >= 0, ...
            for (i, fi) in enumerate(func)
                new_con = @constraint(model, fi >= 0)
                if !isempty(original_name)
                    set_name(new_con, _insert_suffix(original_name, "_$i"))
                end
            end

        elseif set isa MOI.Nonpositives
            # [f1; f2; ...] ∈ Nonpos  →  f1 <= 0, f2 <= 0, ...
            for (i, fi) in enumerate(func)
                new_con = @constraint(model, fi <= 0)
                if !isempty(original_name)
                    set_name(new_con, _insert_suffix(original_name, "_$i"))
                end
            end
        end

        delete(model, con)
    end
end

include("oracleTypicalClassical.jl")
include("oracleTypicalSeparable.jl")
include("oracleTypicalUnified.jl")
include("oracleTypicalPareto.jl")
