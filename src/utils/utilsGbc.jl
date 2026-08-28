"""
    _parse_gbc_result(result, x_vars) -> (gbc_lhs, gbc_rhs, gbc_sense)

Parse and validate the value returned by a user model-update function and
convert it into the internal GBC representation.

# Returns
(gbc_lhs, gbc_rhs, gbc_sense)
- `gbc_lhs :: Vector{VariableRef}`  
  Left-hand-side variables for the generalized bound constraint. Every item
  must be a `VariableRef` that is **not** one of the copied master variables
  supplied in `x_vars`.
- `gbc_rhs :: Vector{Union{VariableRef, AffExpr}}`  
  Right-hand-side expressions. Each item must be either
  - a `VariableRef` drawn from `x_vars` (i.e., a *copied* master variable), or
  - an `AffExpr` expressed only in terms of variables from `x_vars` (and
    numeric constants).
- `gbc_sense :: Vector{GBCBoundType}`  
  A vector of bound senses for each constraint (for example `UpperBound`,
  `LowerBound`, or `Fixed`).

# Arguments
- `result` :: Any
  The raw return value from the user model-update function. This function must
  parse and normalize `result` into the three-tuple described above.
- `x_vars` :: Vector{VariableRef}  
  A vector of submodel variables that are *copies* of the coupling master
  variables. Used only for validation.

# Validation
1. `gbc_lhs` must be a `Vector{VariableRef}`. None of the entries in `gbc_lhs`
   may be any of the `VariableRef`s listed in `x_vars`.
2. `gbc_rhs` must be a `Vector{Union{VariableRef, AffExpr}}`. Every `VariableRef`
   appearing in `gbc_rhs` must be one of the `VariableRef`s in `x_vars`.
3. `gbc_sense` must be a `Vector{GBCBoundType}` with the same length as `gbc_lhs`
   and `gbc_rhs`.

# Errors
If any validation rule is violated, raise a descriptive exception (for example
`ArgumentError`). Error messages should identify:
- which part (`gbc_lhs` / `gbc_rhs` / `gbc_sense`) failed validation, and
- the offending element(s) and why (e.g., "lhs contains copied master variable `x[3]`",
  or "rhs contains variable `y` not in `x_vars`").

# Examples
See the test suite for representative `result` shapes and the exact exception
types.
"""
function _parse_gbc_result(result, x_vars::Vector{VariableRef})
    # No GBC - must explicitly return nothing
    if result === nothing
        return VariableRef[], Union{VariableRef, AffExpr}[], GBCBoundType[]
    end
    
    # Format: (gbc_lhs, gbc_rhs, gbc_sense)
    if length(result) == 3
        gbc_lhs, gbc_rhs, gbc_sense = result
        
        # Validate lengths
        n = length(gbc_lhs)
        if length(gbc_rhs) != n || length(gbc_sense) != n
            throw(DimensionMismatch(
                "All GBC vectors must have the same length. " *
                "Got length(gbc_lhs) = $n, length(gbc_rhs) = $(length(gbc_rhs)), " *
                "length(gbc_sense) = $(length(gbc_sense))."
            ))
        end
        
        if isempty(gbc_lhs)
            @warn "GBC tuple returned but gbc_lhs is empty. No GBC constraints will be applied."
            return VariableRef[], Union{VariableRef, AffExpr}[], GBCBoundType[]
        end
        
        # Build set of valid master variable indices for fast lookup
        x_indices = Set(v.index for v in x_vars)
        
        # Validate gbc_lhs: must be a vector of VariableRef and must NOT contain copied master variables (should be subproblem variables)
        for (i, lhs_var) in enumerate(gbc_lhs)
            if !(lhs_var isa VariableRef)
                throw(ArgumentError("gbc_lhs[$i] is not a VariableRef (got $(typeof(lhs_var))). Each lhs entry must be a VariableRef."))
            end
            if lhs_var.index in x_indices
                throw(ArgumentError(
                    "gbc_lhs[$i] contains a copied master variable $lhs_var. " *
                    "Each item of gbc_lhs should only contain a single subproblem variable. " *
                    "The GBC relation should be: subproblem_var <=/>=/== affine_func(master_vars)."
                ))
            end
        end
        
        # Validate gbc_rhs: must ONLY contain copied master variables in affine form
        for (i, rhs) in enumerate(gbc_rhs)
            if rhs isa VariableRef
                if !(rhs.index in x_indices)
                    throw(ArgumentError(
                        "gbc_rhs[$i] contains a variable $(rhs) that is not a copied master variable. " *
                        "gbc_rhs should only reference copied master variables passed to the model-update function. "
                    ))
                end
            elseif rhs isa AffExpr  # AffExpr
                for (var, _) in rhs.terms
                    if !(var.index in x_indices)
                        throw(ArgumentError(
                            "gbc_rhs[$i] contains an AffExpr with a variable $(var) that is not a copied master variable. " *
                            "gbc_rhs should only reference copied master variables passed to the model-update function. "
                        ))
                    end
                end
            else
                throw(ArgumentError(
                            "gbc_rhs[$i] contains a non-affine expression $rhs. " *
                            "we currently only support affine `gbc_rhs`. "
                        ))
            end
        end

        # Validate gbc_sense: must be a vector of GBCBoundType
        for (i, sense) in enumerate(gbc_sense)
            if !(sense isa GBCBoundType)
                throw(ArgumentError("gbc_sense[$i] is not a GBCBoundType (got $(sense)). Each entry must be among `UpperBound`, `LowerBound`, and `Fixed`."))
            end
        end
        
        return Vector{VariableRef}(gbc_lhs), Vector{Union{VariableRef, AffExpr}}(gbc_rhs), Vector{GBCBoundType}(gbc_sense)
    end
    
    throw(ArgumentError(
        "Invalid return format of model-update function. Expected nothing or (gbc_lhs, gbc_rhs, gbc_sense). " *
        "Got $(length(result))."
    ))
end

"""
    _set_gbc_bounds!(gbc_lhs, gbc_rhs, gbc_sense, x_value)

Set the bounds for GBC LHS variables by evaluating RHS expressions with given x values.
Optimized to handle VariableRef directly without AffExpr overhead.
"""
function _set_gbc_bounds!(gbc_lhs::Vector{VariableRef}, 
                          gbc_rhs::Vector{Union{VariableRef, AffExpr}},
                          gbc_sense::Vector{GBCBoundType},
                          x_value::Vector{Float64})
    for i in 1:length(gbc_lhs)
        rhs = gbc_rhs[i]
        # Evaluate RHS expression
        if rhs isa VariableRef
            # Fast path for simple variable bounds
            bound_value = x_value[rhs.index.value]
        else
            # General path for affine expressions
            bound_value = value(v -> x_value[v.index.value], rhs)
        end
        
        # Set bound based on sense
        if gbc_sense[i] == UpperBound
            set_upper_bound(gbc_lhs[i], bound_value)
        elseif gbc_sense[i] == LowerBound
            set_lower_bound(gbc_lhs[i], bound_value)
        else  # Fixed
            fix(gbc_lhs[i], bound_value; force=true)
        end
    end
end

"""
    _accumulate_gbc_duals!(a_x, gbc_lhs, gbc_rhs, gbc_sense)

Accumulate dual values from GBC constraints into the cut coefficients.
Optimized to handle VariableRef directly without AffExpr overhead.
"""
function _accumulate_gbc_duals!(a_x::Vector{Float64},
                                gbc_lhs::Vector{VariableRef},
                                gbc_rhs::Vector{Union{VariableRef, AffExpr}},
                                gbc_sense::Vector{GBCBoundType})
    for i in 1:length(gbc_lhs)
        # Get dual value based on bound sense
        if gbc_sense[i] == UpperBound
            dual_val = dual(UpperBoundRef(gbc_lhs[i]))
        elseif gbc_sense[i] == LowerBound
            dual_val = dual(LowerBoundRef(gbc_lhs[i]))
        else  # Fixed
            dual_val = dual(FixRef(gbc_lhs[i]))
        end
        
        # Accumulate to corresponding x positions
        rhs = gbc_rhs[i]
        if rhs isa VariableRef
            # Fast path for simple variable bounds
            a_x[rhs.index.value] += dual_val
        else
            # General path for affine expressions
            for (var, coef) in rhs.terms
                a_x[var.index.value] += coef * dual_val
            end
        end
    end
end