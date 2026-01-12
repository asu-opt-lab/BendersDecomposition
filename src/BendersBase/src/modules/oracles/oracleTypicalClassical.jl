export ClassicalOracle, ClassicalOracleParam

const ClassicalOracleParam = BasicOracleParam

mutable struct ClassicalOracle <: AbstractTypicalOracle
    
    param::ClassicalOracleParam

    model::Model
    fixed_x_constraints::Vector{ConstraintRef}

    gbc_y::Vector{VariableRef}
    gbc_x_indices::Vector{Vector{Int}}
    gbc_coefficients::Vector{Vector{Float64}}
    gbc_bound_type::Vector{GBCBoundType}


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

            # Build mapping from variable index to x vector position
            idx_to_pos = Dict{Int,Int}()
            for (pos, v) in enumerate(x)
                vi = JuMP.index(v)
                idx_to_pos[vi.value] = pos
            end

            # Build the submodel using user-defined customization, passing the copied variables
            result = customize(model, data, scen_idx; x_copy...)
            
            # Parse the result to extract GBC information
            gbc_y, gbc_x_indices, gbc_coefficients, gbc_bound_type = _parse_gbc_result(result, idx_to_pos)

            new(param, model, fix_x, gbc_y, gbc_x_indices, gbc_coefficients, gbc_bound_type)
    end

    ClassicalOracle() = new()
end

"""
    _parse_gbc_result(result, idx_to_pos) -> (gbc_y, gbc_x_indices, gbc_coefficients, gbc_bound_type)

Parse the result returned by customize function and convert to internal GBC representation.

Supported formats:
- `nothing` or no tuple: No GBC constraints
- `(gbc_y, gbc_x)`: Legacy format, one-to-one mapping with coefficient 1.0, default UpperBound
- `(gbc_y, gbc_x_indices, gbc_coefficients, gbc_bound_type)`: New format with full specification
"""
function _parse_gbc_result(result, idx_to_pos::Dict{Int,Int})
    # No GBC
    if result === nothing || !(result isa Tuple)
        return VariableRef[], Vector{Int}[], Vector{Float64}[], GBCBoundType[]
    end
    
    # Legacy format: (gbc_y, gbc_x)
    if length(result) == 2
        gbc_y, gbc_x = result
        if length(gbc_y) != length(gbc_x)
            throw(DimensionMismatch(
                "gbc_y and gbc_x returned by customize function must have the same length. " *
                "Got length(gbc_y) = $(length(gbc_y)), length(gbc_x) = $(length(gbc_x)). " *
                "Each gbc_y[i] should correspond to gbc_x[i] for the Generalized Bound Constraints."
            ))
        end
        
        if isempty(gbc_y)
            return VariableRef[], Vector{Int}[], Vector{Float64}[], GBCBoundType[]
        end
        
        # Convert to new format: one-to-one mapping with coefficient 1.0, default UpperBound
        gbc_x_indices = [[idx_to_pos[JuMP.index(v).value]] for v in gbc_x]
        gbc_coefficients = [[1.0] for _ in gbc_x]
        gbc_bound_type = fill(UpperBound, length(gbc_y))
        
        return gbc_y, gbc_x_indices, gbc_coefficients, gbc_bound_type
    end
    
    # New format: (gbc_y, gbc_x_indices, gbc_coefficients, gbc_bound_type)
    if length(result) == 4
        gbc_y, gbc_x_indices_raw, gbc_coefficients, gbc_bound_type = result
        
        # Validate lengths
        n = length(gbc_y)
        if length(gbc_x_indices_raw) != n || length(gbc_coefficients) != n || length(gbc_bound_type) != n
            throw(DimensionMismatch(
                "All GBC vectors must have the same length. " *
                "Got length(gbc_y) = $n, length(gbc_x_indices) = $(length(gbc_x_indices_raw)), " *
                "length(gbc_coefficients) = $(length(gbc_coefficients)), length(gbc_bound_type) = $(length(gbc_bound_type))."
            ))
        end
        
        # Validate inner lengths match
        for i in 1:n
            if length(gbc_x_indices_raw[i]) != length(gbc_coefficients[i])
                throw(DimensionMismatch(
                    "For GBC entry $i, gbc_x_indices and gbc_coefficients must have the same length. " *
                    "Got $(length(gbc_x_indices_raw[i])) and $(length(gbc_coefficients[i]))."
                ))
            end
        end
        
        if isempty(gbc_y)
            return VariableRef[], Vector{Int}[], Vector{Float64}[], GBCBoundType[]
        end
        
        # Convert variable references to indices if needed
        gbc_x_indices = Vector{Vector{Int}}(undef, n)
        for i in 1:n
            if eltype(gbc_x_indices_raw[i]) <: VariableRef
                gbc_x_indices[i] = [idx_to_pos[JuMP.index(v).value] for v in gbc_x_indices_raw[i]]
            else
                gbc_x_indices[i] = collect(Int, gbc_x_indices_raw[i])
            end
        end
        
        return gbc_y, gbc_x_indices, Vector{Vector{Float64}}(gbc_coefficients), Vector{GBCBoundType}(gbc_bound_type)
    end
    
    throw(ArgumentError(
        "Invalid GBC result format. Expected nothing, (gbc_y, gbc_x), or " *
        "(gbc_y, gbc_x_indices, gbc_coefficients, gbc_bound_type). Got tuple of length $(length(result))."
    ))
end

function generate_cuts(oracle::ClassicalOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    set_time_limit_sec(oracle.model, time_limit)
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    
    # Set GBC bounds based on affine combination and bound type
    _set_gbc_bounds!(oracle.gbc_y, oracle.gbc_x_indices, oracle.gbc_coefficients, 
                     oracle.gbc_bound_type, x_value)
    
    optimize!(oracle.model)
    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end
    
    status = dual_status(oracle.model)
    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(oracle.model)

        a_x = dual.(oracle.fixed_x_constraints)
        
        # Accumulate GBC dual values
        _accumulate_gbc_duals!(a_x, oracle.gbc_y, oracle.gbc_x_indices, 
                              oracle.gbc_coefficients, oracle.gbc_bound_type)
        
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
            _accumulate_gbc_duals!(a_x, oracle.gbc_y, oracle.gbc_x_indices, 
                                  oracle.gbc_coefficients, oracle.gbc_bound_type)
            
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
    _set_gbc_bounds!(gbc_y, gbc_x_indices, gbc_coefficients, gbc_bound_type, x_value)

Set the bounds for GBC variables based on affine combination of x values and bound type.
"""
function _set_gbc_bounds!(gbc_y::Vector{VariableRef}, 
                          gbc_x_indices::Vector{Vector{Int}},
                          gbc_coefficients::Vector{Vector{Float64}},
                          gbc_bound_type::Vector{GBCBoundType},
                          x_value::Vector{Float64})
    for i in 1:length(gbc_y)
        # Calculate affine combination: sum(coeff[j] * x_value[idx[j]])
        bound_value = sum(gbc_coefficients[i][j] * x_value[gbc_x_indices[i][j]] 
                         for j in 1:length(gbc_x_indices[i]))
        
        # Set bound based on type
        if gbc_bound_type[i] == UpperBound
            set_upper_bound(gbc_y[i], bound_value)
        elseif gbc_bound_type[i] == LowerBound
            set_lower_bound(gbc_y[i], bound_value)
        else  # FixedBound
            fix(gbc_y[i], bound_value; force=true)
        end
    end
end

"""
    _accumulate_gbc_duals!(a_x, gbc_y, gbc_x_indices, gbc_coefficients, gbc_bound_type)

Accumulate dual values from GBC constraints into the cut coefficients.
For each GBC entry: a_x[idx[j]] += coeff[j] * dual(BoundRef(y))
"""
function _accumulate_gbc_duals!(a_x::Vector{Float64},
                                gbc_y::Vector{VariableRef},
                                gbc_x_indices::Vector{Vector{Int}},
                                gbc_coefficients::Vector{Vector{Float64}},
                                gbc_bound_type::Vector{GBCBoundType})
    for i in 1:length(gbc_y)
        # Get dual value based on bound type
        if gbc_bound_type[i] == UpperBound
            dual_val = dual(UpperBoundRef(gbc_y[i]))
        elseif gbc_bound_type[i] == LowerBound
            dual_val = dual(LowerBoundRef(gbc_y[i]))
        else  # FixedBound
            dual_val = dual(FixRef(gbc_y[i]))
        end
        
        # Accumulate to corresponding x positions with coefficients
        for j in 1:length(gbc_x_indices[i])
            a_x[gbc_x_indices[i][j]] += gbc_coefficients[i][j] * dual_val
        end
    end
end




