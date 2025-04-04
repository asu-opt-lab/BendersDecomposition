export GenericSequential
export solve!

struct GenericSequential <: AbstractSequential
    time_limit::Float64
    iteration_limit::Int
    gap_tolerance::Float64
    verbose::Bool
end


function solve!(env::AbstractBendersEnv, SolutionProcedure::GenericSequential, cut_strategy::AbstractCutStrategy)
    # parameters
    time_limit = SolutionProcedure.time_limit
    iteration_limit = SolutionProcedure.iteration_limit
    gap_tolerance = SolutionProcedure.gap_tolerance
    verbose = SolutionProcedure.verbose

    state = BendersState()
    start_time = time()
    while true
        state.iteration += 1
        
        # Solve master problem
        master_time = @elapsed begin
            solve_master!(env.master)
            state.LB = env.master.objective_value
        end
        state.master_time += master_time
        
        # Solve sub problem
        sub_time = @elapsed begin
            solve_sub!(env.sub, env.master.integer_variable_values)
            cuts, sub_obj_val = generate_cuts(env, cut_strategy)
            update_upper_bound_and_gap!(state, env, sub_obj_val)
        end
        state.sub_time += sub_time
        state.total_time = time() - start_time

        # Update state and record information
        verbose && print_iteration_info(state)

        # Check termination criteria
        is_terminated(state, iteration_limit, time_limit, gap_tolerance) && break

        # Generate and add cuts
        add_cuts!(env, cuts)
    end
    
end




struct GenericCallback <: AbstractCallback 
    time_limit::Union{Float64,Nothing}
    gap_tolerance::Union{Float64}
    lazy_callback::Function
    user_callback::Union{Function, Nothing}
    verbose::Bool
end

# Solve with GenericCallback
function solve!(env::AbstractBendersEnv, solution_procedure::GenericCallback, cut_strategy::AbstractCutStrategy)
    # Set up lazy callback
    function lazy_callback_wrapper(cb_data)
        solution_procedure.lazy_callback(cb_data, env, cut_strategy)
    end
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback_wrapper)
    
    # Set up user callback if provided
    if solution_procedure.user_callback !== nothing
        function user_callback_wrapper(cb_data)
            solution_procedure.user_callback(cb_data, env, cut_strategy)
        end
        set_attribute(env.master.model, MOI.UserCutCallback(), user_callback_wrapper)
    end

    # Configure solver parameters
    set_time_limit_sec(env.master.model, solution_procedure.time_limit)
    set_optimizer_attribute(env.master.model, MOI.Silent(), false)

    # Solve and log results
    start_time = time()
    JuMP.optimize!(env.master.model)
    
    # Log performance metrics
    @info "Optimization results:" begin
        "node count" => JuMP.node_count(env.master.model)
        "elapsed time" => time() - start_time
        "objective bound" => JuMP.objective_bound(env.master.model)
        "objective value" => JuMP.objective_value(env.master.model)
        "relative gap" => JuMP.relative_gap(env.master.model)
    end
end

function lazy_callback(cb_data, env::AbstractBendersEnv, cut_strategy::AbstractCutStrategy)
    # Only add cuts at integer nodes
    status = JuMP.callback_node_status(cb_data, env.master.model)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER    
        # Get current solution values
        env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
        env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])
        
        # Generate cuts based on current solution
        solve_sub!(env.sub, env.master.x_value)
        cuts, sub_obj_value = generate_cuts(env, cut_strategy)
        
        # Add violated cuts
        add_cuts!(env, cuts, sub_obj_value, cb_data)
    end
end

# Base add_cuts! function
function add_cuts!(env::AbstractBendersEnv, expressions::Union{Vector{Any}, Any}, sub_obj_values::Union{Vector{Float64}, Float64}, cb_data)
    throw(ArgumentError("Unsupported types for add_cuts!: expressions=$(typeof(expressions)), sub_obj_values=$(typeof(sub_obj_values))"))
end

# Add cuts for vector inputs
function add_cuts!(env::AbstractBendersEnv, expressions::Vector{Any}, sub_obj_values::Vector{Float64}, cb_data)
    for (idx, (expr, sub_obj)) in enumerate(zip(expressions, sub_obj_values))
        # Only add cut if it's violated (with numerical tolerance)
        if env.master.t_value[idx] <= sub_obj - 1e-06
            cut = @build_constraint(0 >= expr)
            MOI.submit(env.master.model, MOI.LazyConstraint(cb_data), cut)
        end
    end
end

# Add cuts for single expression
function add_cuts!(env::AbstractBendersEnv, expr::Any, sub_obj_value::Float64, cb_data)
    # Only add cut if it's violated (with numerical tolerance)
    if env.master.t_value <= sub_obj_value - 1e-06
        cut = @build_constraint(0 >= expr)
        MOI.submit(env.master.model, MOI.LazyConstraint(cb_data), cut)
    end
end

