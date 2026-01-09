export ClassicalOracle, ClassicalOracleParam

const ClassicalOracleParam = BasicOracleParam

mutable struct ClassicalOracle <: AbstractTypicalOracle
    
    param::ClassicalOracleParam

    model::Model
    fixed_x_constraints::Vector{ConstraintRef}

    gbc_y::Vector{VariableRef}
    gbc_x_idx::Vector{Int} 


    function ClassicalOracle(data::AbstractData, master::Master; 
                            customize = customize_sub_model!,
                            scen_idx::Int=0, 
                            param::ClassicalOracleParam = ClassicalOracleParam())
    
            @debug "Building classical oracle"
            model = Model()

            # Copy the master’s coupling variables into the submodel (with identical axes and symbols)
            x_copy = copy_variables!(model, master.x_tuple)

            # Build the submodel using user-defined customization, passing the copied variables
            result = customize(model, data, scen_idx; x_copy...)
            if result isa Tuple && length(result) == 2
                gbc_y, gbc_x = result
            else
                gbc_y = VariableRef[]
                gbc_x = VariableRef[]
            end

            # Collect all copied master variables and add linking constraint
            x = var_from_tuple(x_copy)
            @constraint(model, fix_x, x .== 0)

            # Build mapping from gbc_x[i] to x vector position
            if !isempty(gbc_x)
                idx_to_pos = Dict{Int,Int}()
                for (pos, v) in enumerate(x)
                    vi = JuMP.index(v)
                    idx_to_pos[vi.value] = pos
                end
                gbc_x_idx = Int[idx_to_pos[JuMP.index(v).value] for v in gbc_x]
            else
                gbc_x_idx = Int[]
            end
            new(param, model, fix_x, gbc_y, gbc_x_idx)
    end

    ClassicalOracle() = new()
end

function generate_cuts(oracle::ClassicalOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    set_time_limit_sec(oracle.model, time_limit)
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    if !isempty(oracle.gbc_y)
        for i in 1:length(oracle.gbc_y)
            set_upper_bound(oracle.gbc_y[i], x_value[oracle.gbc_x_idx[i]])
        end
    end
    optimize!(oracle.model)
    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end
    
    status = dual_status(oracle.model)
    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(oracle.model)

        a_x = dual.(oracle.fixed_x_constraints)
        if !isempty(oracle.gbc_y)
            for i in 1:length(oracle.gbc_y)
                a_x[oracle.gbc_x_idx[i]] += dual(UpperBoundRef(oracle.gbc_y[i]))
            end
        end
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
            if !isempty(oracle.gbc_y)
                for i in 1:length(oracle.gbc_y)
                    a_x[oracle.gbc_x_idx[i]] += dual(UpperBoundRef(oracle.gbc_y[i]))
                end
            end
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




