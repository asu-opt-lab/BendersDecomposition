export ClassicalOracle, ClassicalOracleParam

const ClassicalOracleParam = BasicOracleParam

mutable struct ClassicalOracle <: AbstractTypicalOracle
    
    param::ClassicalOracleParam

    model::Model
    fixed_x_constraints::Vector{ConstraintRef}

    gbc_lhs::Vector{VariableRef}
    gbc_rhs::Vector{Union{VariableRef, AffExpr}}
    gbc_sense::Vector{GBCBoundType}


    function ClassicalOracle(data::AbstractData, master::Master; 
                            customize = customize_sub_model!,
                            scen_idx::Int=0, 
                            param::ClassicalOracleParam = ClassicalOracleParam())
    
            @debug "Building classical oracle"
            model = Model()

            # Copy the master's coupling variables into the submodel (with identical axes and symbols)
            x_copy = copy_variables!(model, master.x_tuple)

            # Collect all copied master variables and add linking constraint
            x = var_from_tuple(x_copy)
            @constraint(model, fix_x, x .== 0)

            # Build the submodel using user-defined customization, passing the copied variables
            result = customize(model, data, scen_idx; x_copy...)
            
            # Parse the result to extract GBC information
            gbc_lhs, gbc_rhs, gbc_sense = _parse_gbc_result(result, x)

            new(param, model, fix_x, gbc_lhs, gbc_rhs, gbc_sense)
    end

    ClassicalOracle() = new()
end

"""
    _parse_gbc_result(result, x_vars) -> (gbc_lhs, gbc_rhs, gbc_sense)

Parse the result returned by customize function and convert to internal GBC representation.
Returns LHS variables, RHS expressions (AffExpr), and bound sense.

Supported formats:
- `nothing` or no tuple: No GBC constraints
- `(gbc_lhs, gbc_rhs, gbc_sense)`: gbc_rhs is Vector{Union{VariableRef, AffExpr}}

# Arguments
- `result`: The return value from the customize function
- `x_vars`: Vector of master variables (copied into submodel) for validation

# Validation
- gbc_lhs must contain subproblem variables (NOT master variables from x)
- gbc_rhs must contain only master variables (from x)
"""
function _parse_gbc_result(result, x_vars::Vector{VariableRef})
    # No GBC
    if result === nothing || !(result isa Tuple)
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
        
        # Validate gbc_lhs: must NOT contain master variables (should be subproblem variables)
        for (i, lhs_var) in enumerate(gbc_lhs)
            if lhs_var.index in x_indices
                throw(ArgumentError(
                    "gbc_lhs[$i] contains a master variable. " *
                    "gbc_lhs should only contain subproblem variables. " *
                    "The GBC relation should be: subproblem_var <=/>=/== f(master_vars)."
                ))
            end
        end
        
        # Validate gbc_rhs: must ONLY contain master variables
        for (i, rhs) in enumerate(gbc_rhs)
            if rhs isa VariableRef
                if !(rhs.index in x_indices)
                    throw(ArgumentError(
                        "gbc_rhs[$i] contains a variable that is not a registered master variable. " *
                        "gbc_rhs should only reference master variables passed to customize(). " *
                        "Got variable with index $(rhs.index.value), but valid master indices are: $(sort([v.value for v in x_indices]))."
                    ))
                end
            else  # AffExpr
                for (var, _) in rhs.terms
                    if !(var.index in x_indices)
                        throw(ArgumentError(
                            "gbc_rhs[$i] contains an AffExpr with a variable that is not a registered master variable (x). " *
                            "gbc_rhs should only reference master variables passed to customize(). " *
                            "Got variable with index $(var.index.value), but valid master indices are: $(sort([v.value for v in x_indices]))."
                        ))
                    end
                end
            end
        end
        
        return Vector{VariableRef}(gbc_lhs), Vector{Union{VariableRef, AffExpr}}(gbc_rhs), Vector{GBCBoundType}(gbc_sense)
    end
    
    throw(ArgumentError(
        "Invalid GBC result format. Expected nothing or (gbc_lhs, gbc_rhs, gbc_sense). " *
        "Got tuple of length $(length(result))."
    ))
end


function generate_cuts(oracle::ClassicalOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    set_time_limit_sec(oracle.model, time_limit)
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    
    # Set GBC bounds based on expression evaluation
    _set_gbc_bounds!(oracle.gbc_lhs, oracle.gbc_rhs, oracle.gbc_sense, x_value)
    
    optimize!(oracle.model)
    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end
    
    status = dual_status(oracle.model)
    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(oracle.model)

        a_x = dual.(oracle.fixed_x_constraints)
        
        # Accumulate GBC dual values
        _accumulate_gbc_duals!(a_x, oracle.gbc_lhs, oracle.gbc_rhs, oracle.gbc_sense)
        
        a_t = [-1.0] 
        a_0 = sub_obj_val - a_x'*x_value 
        if sub_obj_val >= t_value[1] * (1 + oracle.param.rtol) + oracle.param.atol / tol_normalize
            return false, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        else
            return true, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        end

    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(oracle.model)
            dual_sub_obj_val = dual_objective_value(oracle.model)
            @debug "dual_sub_obj_val = $dual_sub_obj_val"
            a_x = dual.(oracle.fixed_x_constraints)
            
            # Accumulate GBC dual values
            _accumulate_gbc_duals!(a_x, oracle.gbc_lhs, oracle.gbc_rhs, oracle.gbc_sense)
            
            a_t = [0.0]
            a_0 = dual_sub_obj_val - a_x' * x_value 
            if dual_sub_obj_val >= oracle.param.zero_tol / tol_normalize
                return false, [Hyperplane(a_x, a_t, a_0)], [Inf]
            else
                return true, [Hyperplane(a_x, a_t, a_0)], [Inf]
            end
        end
    else
        throw(UnexpectedModelStatusException("ClassicalOracle: $(status). This is likely a numerical issue. Please try using other oracles, such as unified oracle or pareto oracle."))
    end
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
        else  # FixedBound
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
        else  # FixedBound
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




