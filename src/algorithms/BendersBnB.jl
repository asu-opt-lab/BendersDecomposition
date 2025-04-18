export BendersBnB, solve!

include("callback.jl")

mutable struct BendersBnB <: AbstractBendersCallback
    data::Data
    master::AbstractMaster 

    param::BendersBnBParam 

    root_preprocessing::Union{RootNodePreprocessing, Nothing} 
    lazy_callback::AbstractLazyCallback
    user_callback::Union{AbstractUserCallback, Nothing} 

    obj_value::Float64 
    termination_status::TerminationStatus 

    function BendersBnB(data; param::BendersBnBParam = BendersBnBParam())
        new(data, Master(data), param, RootNodePreprocessing(data), LazyCallback(data), UserCallback(data), Inf, NotSolved())
    end

    function BendersBnB(data, master::AbstractMaster, root_preprocessing::Union{RootNodePreprocessing, Nothing}, lazy_callback::AbstractLazyCallback, user_callback::Union{AbstractUserCallback, Nothing}; param::BendersBnBParam = BendersBnBParam())
        new(data, master, param, root_preprocessing, lazy_callback, user_callback, Inf, NotSolved())
    end
end

function solve!(env::BendersBnB) 
    log = BendersBnBLog()
    param = env.param
    start_time = time()
    
    if env.root_preprocessing !== nothing
        root_node_time = root_node_processing!(env.data, env.master, env.root_preprocessing)
    end
    
    function lazy_callback_wrapper(cb_data)
        lazy_callback(cb_data, env.master.model, log, env.lazy_callback)
    end
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback_wrapper)
    
    if env.user_callback !== nothing
        function user_callback_wrapper(cb_data)
            user_callback(cb_data, env.master.model, log, env.user_callback)
        end
        set_attribute(env.master.model, MOI.UserCutCallback(), user_callback_wrapper)
    end
    
    set_time_limit_sec(env.master.model, param.time_limit)
    set_optimizer_attribute(env.master.model, MOI.Silent(), !param.verbose)
    set_optimizer_attribute(env.master.model, MOI.RelativeGapTolerance(), param.gap_tolerance)
    
    JuMP.optimize!(env.master.model)
    
    status = termination_status(env.master.model)
    if status == MOI.OPTIMAL
        env.termination_status = Optimal()
        env.obj_value = JuMP.objective_value(env.master.model)
    elseif status == MOI.TIME_LIMIT
        env.termination_status = TimeLimit()
        env.obj_value = has_values(env.master.model) ? JuMP.objective_value(env.master.model) : Inf
    else
        env.termination_status = InfeasibleOrNumericalIssue()
        env.obj_value = Inf
    end
    
    elapsed_time = time() - start_time
    
    if param.verbose
        @info "Node count: $(JuMP.node_count(env.master.model))"
        @info "Elapsed time: $(elapsed_time)"
        @info "Objective bound: $(JuMP.objective_bound(env.master.model))"
        @info "Objective value: $(env.obj_value)"
        @info "Relative gap: $(JuMP.relative_gap(env.master.model))"
    end
    
    return env.obj_value, elapsed_time
end



