export BendersCallback, solve!

mutable struct BendersCallback <: AbstractBendersCallback
    data::Data
    master::AbstractMaster
    oracle::AbstractOracle

    param::BendersCallbackParam # parameters for the algorithm
    lazy_callback::Function
    user_callback::Union{Function, Nothing}
    # result
    obj_value::Float64
    termination_status::TerminationStatus

    function BendersCallback(data, master::AbstractMaster, oracle::AbstractOracle, lazy_callback::Function, user_callback::Union{Function, Nothing}; param::BendersCallbackParam = BendersCallbackParam()) 
        
        # Constructor with master and oracle provided
        env = new(data, master, oracle, param, lazy_callback, user_callback, Inf, NotSolved())

    end

    function BendersCallback(data, master::AbstractMaster, oracle::AbstractOracle; param::BendersCallbackParam = BendersCallbackParam()) 
        
        # Constructor with master and oracle provided
        env = new(data, master, oracle, param, default_lazy_callback, nothing, Inf, NotSolved())

    end

    function BendersCallback(data; param::BendersCallbackParam = BendersCallbackParam())
        # Constructor with just data, creating default master and oracle
        new(data, Master(data), ClassicalOracle(data), param, default_lazy_callback, nothing, Inf, NotSolved())
    end
end

"""
Run BendersCallback decomposition using lazy constraint callbacks
"""
function solve!(env::BendersCallback) 
    param = env.param
    start_time = time()
    
    if param.preprocessing_type !== nothing
        root_node_time = root_node_processing!(env, param.preprocessing_type)
    end
    # Pre-processing step based on the type of preprocessing selected
    # if !(param.preprocessing_type isa NoPreprocessing)
    #     # Relax integrality for root node processing
    #     relaxed_vars = relax_integrality(env.master.model)
        
    #     # Run preprocessing at root node to add initial cuts
    #     root_node_time = @elapsed begin
    #         # Get a copy of the root node parameters with adjusted time limit
    #         root_param = deepcopy(param.root_param)
            
    #         # Run appropriate preprocessing algorithm based on type
    #         if param.preprocessing_type isa SeqPreprocessing
    #             # Run sequential Benders at root node
    #             seq_solver = BendersSeq(env.data, env.master, env.oracle; param=root_param)
    #             solve!(seq_solver)
    #         elseif param.preprocessing_type isa SeqInOutPreprocessing
    #             # Run stabilized in-out Benders at root node
    #             # Use default starting point (or could be customized)
    #             inout_solver = BendersSeqInOut(env.data, env.master, env.oracle; param=root_param)
    #             solve!(inout_solver)
    #         end

    #     end
        
    #     # Restore integrality
    #     set_binary.(env.master.model[:x])
        
    #     if param.verbose
    #         @info "Root node preprocessing completed in $(root_node_time) seconds"
    #         @info "Remaining time limit: $(param.time_limit) seconds"
    #     end
    # end
    
    # Define the lazy constraint callback function
    
    function lazy_callback_wrapper(cb_data)
        env.lazy_callback(cb_data, env)
    end
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback_wrapper)

    # function lazy_callback(cb_data)
    #     status = JuMP.callback_node_status(cb_data, env.master.model)
    #     if status == MOI.CALLBACK_NODE_STATUS_INTEGER
    #         # Get current values from the master problem
    #         x_vals = JuMP.callback_value.(cb_data, env.master.model[:x])
    #         t_vals = JuMP.callback_value.(cb_data, env.master.model[:t])

    #         # Generate cuts based on current solution
    #         is_in_L, hyperplanes, f_x = generate_cuts(env.oracle, x_vals, t_vals)
            
    #         # If not feasible, add cuts
    #         cuts = !is_in_L ? hyperplanes_to_expression(env.master.model, hyperplanes, env.master.model[:x], env.master.model[:t]) : []

    #         if !isempty(cuts)
    #             for cut in cuts
    #                 cut_constraint = @build_constraint(0 >= cut)
    #                 MOI.submit(env.master.model, MOI.LazyConstraint(cb_data), cut_constraint)
    #             end
    #         end
    #     end
    # end
    # set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback)

    
    # Set up user callback if provided
    if env.user_callback !== nothing
        function user_callback_wrapper(cb_data)
            env.user_callback(cb_data, env)
        end
        set_attribute(env.master.model, MOI.UserCutCallback(), user_callback_wrapper)
    end
    
    # Set solver parameters
    set_time_limit_sec(env.master.model, param.time_limit)
    set_optimizer_attribute(env.master.model, MOI.Silent(), !param.verbose)
    set_optimizer_attribute(env.master.model, MOI.RelativeGapTolerance(), param.gap_tolerance)

    
    # Solve the master problem
    JuMP.optimize!(env.master.model)
    
    # Record results
    if termination_status(env.master.model) == MOI.OPTIMAL
        env.termination_status = Optimal()
        env.obj_value = JuMP.objective_value(env.master.model)
    elseif termination_status(env.master.model) == MOI.TIME_LIMIT
        env.termination_status = TimeLimit()
        env.obj_value = has_values(env.master.model) ? JuMP.objective_value(env.master.model) : Inf
    else
        env.termination_status = InfeasibleOrNumericalIssue()
        env.obj_value = Inf
    end
    
    # Print summary if verbose
    if param.verbose
        @info "Node count: $(JuMP.node_count(env.master.model))"
        @info "Elapsed time: $(time() - start_time)"
        @info "Objective bound: $(JuMP.objective_bound(env.master.model))"
        @info "Objective value: $(env.obj_value)"
        @info "Relative gap: $(JuMP.relative_gap(env.master.model))"
    end
    
    return env.obj_value, time() - start_time
end 
    
    
function default_lazy_callback(cb_data, env::BendersCallback)
    status = JuMP.callback_node_status(cb_data, env.master.model)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER
        # Get current values from the master problem
        x_vals = JuMP.callback_value.(cb_data, env.master.model[:x])
        t_vals = JuMP.callback_value.(cb_data, env.master.model[:t])

        # Generate cuts based on current solution
        is_in_L, hyperplanes, f_x = generate_cuts(env.oracle, x_vals, t_vals)
        
        # If not feasible, add cuts
        cuts = !is_in_L ? hyperplanes_to_expression(env.master.model, hyperplanes, env.master.model[:x], env.master.model[:t]) : []

        if !isempty(cuts)
            for cut in cuts
                cut_constraint = @build_constraint(0 >= cut)
                MOI.submit(env.master.model, MOI.LazyConstraint(cb_data), cut_constraint)
            end
        end
    end
end

function root_node_processing!(env::BendersCallback, BendersRootSeqType::Type{T}) where T <: AbstractBendersSeq
    
    root_param = deepcopy(env.param.root_param)

    # Relax integrality
    relax_integrality(env.master.model)
    
    # Run preprocessing at root node to add initial cuts
    root_node_time = @elapsed begin
        BendersRootSeq = BendersRootSeqType(env.data, env.master, env.oracle; param=root_param)
        solve!(BendersRootSeq)
    end
    set_binary.(env.master.model[:x])
    return root_node_time
end
