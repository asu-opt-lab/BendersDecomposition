export BendersBnB, solve!
export LazyCallbackParam, UserCallbackParam
export default_lazy_callback, default_user_callback

Base.@kwdef struct LazyCallbackParam
    func::Function = default_lazy_callback
    params::Dict{String, Any} = Dict{String, Any}()
end


Base.@kwdef struct UserCallbackParam
    func::Function = default_user_callback
    params::Dict{String, Any} = Dict{String, Any}("frequency" => 50)
end

Base.@kwdef mutable struct BendersBnB <: AbstractBendersCallback
    data::Data
    master::AbstractMaster = Master(data)
    oracle::AbstractOracle = ClassicalOracle(data)

    param::BendersBnBParam = BendersBnBParam()
    lazy_callback::LazyCallbackParam = LazyCallbackParam()
    user_callback::Union{UserCallbackParam, Nothing} = nothing

    obj_value::Float64 = Inf
    termination_status::TerminationStatus = NotSolved()
end

function solve!(env::BendersBnB) 
    log = BendersBnBLog()
    param = env.param
    start_time = time()
    
    if param.preprocessing_type !== nothing
        root_node_time = root_node_processing!(env, param.preprocessing_type)
    end
    
    function lazy_callback_wrapper(cb_data)
        env.lazy_callback.func(cb_data, env, log)
    end
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback_wrapper)
    
    if env.user_callback !== nothing
        function user_callback_wrapper(cb_data)
            env.user_callback.func(cb_data, env, log)
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

function root_node_processing!(env::BendersBnB, BendersRootSeqType::Type{T}) where T <: AbstractBendersSeq
    root_param = deepcopy(env.param.root_param)

    undo = relax_integrality(env.master.model)
    
    root_node_time = @elapsed begin
        BendersRootSeq = BendersRootSeqType(env.data, env.master, env.oracle; param=root_param)
        solve!(BendersRootSeq)
    end
    
    undo()
    
    return root_node_time
end

function default_lazy_callback(cb_data, env::BendersBnB, log::BendersBnBLog)
    status = JuMP.callback_node_status(cb_data, env.master.model)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER

        state = BendersBnBState()
        if solver_name(env.master.model) == "CPLEX"
            n_count = Ref{CPXINT}()
            ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)
            state.node = n_count[]
        end

        state.values[:x] = JuMP.callback_value.(cb_data, env.master.model[:x])
        state.values[:t] = JuMP.callback_value.(cb_data, env.master.model[:t])

        state.oracle_time = @elapsed begin
            state.is_in_L, hyperplanes, state.f_x = generate_cuts(env.oracle, state.values[:x], state.values[:t])
            cuts = !state.is_in_L ? hyperplanes_to_expression(env.master.model, hyperplanes, env.master.model[:x], env.master.model[:t]) : []
        end

        if !isempty(cuts)
            for cut in cuts
                cut_constraint = @build_constraint(0 >= cut)
                MOI.submit(env.master.model, MOI.LazyConstraint(cb_data), cut_constraint)
                state.num_cuts += 1
            end
        end
        record_node!(log, state, true)
    end
end


function default_user_callback(cb_data, env::BendersBnB, log::BendersBnBLog)
    status = JuMP.callback_node_status(cb_data, env.master.model)
    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        log.num_of_fraction_node += 1
        if log.num_of_fraction_node >= env.user_callback.params["frequency"]
            log.num_of_fraction_node = 0
            state = BendersBnBState()
            if solver_name(env.master.model) == "CPLEX"
                n_count = Ref{CPXINT}()
                ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)
                state.node = n_count[]
            end
            state.values[:x] = JuMP.callback_value.(cb_data, env.master.model[:x])
            state.values[:t] = JuMP.callback_value.(cb_data, env.master.model[:t])
            state.oracle_time = @elapsed begin
                state.is_in_L, hyperplanes, state.f_x = generate_cuts(env.oracle, state.values[:x], state.values[:t])
                cuts = !state.is_in_L ? hyperplanes_to_expression(env.master.model, hyperplanes, env.master.model[:x], env.master.model[:t]) : []
            end
    
            if !isempty(cuts)
                for cut in cuts
                    cut_constraint = @build_constraint(0 >= cut)
                    MOI.submit(env.master.model, MOI.UserCut(cb_data), cut_constraint)
                    state.num_cuts += 1
                end
            end
            record_node!(log, state, false)
        end
    end
end