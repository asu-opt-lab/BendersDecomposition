export RootNodePreprocessing, LazyCallback, UserCallback
export lazy_callback, user_callback
export EmptyCallbackParam, UserCallbackParam


abstract type AbstractCallbackParam end


struct RootNodePreprocessing 
    oracle::AbstractOracle
    seq_type::Type{<:AbstractBendersSeq}
    params::BendersSeqParam

    function RootNodePreprocessing(oracle::AbstractOracle, seq_type::Type{<:AbstractBendersSeq}, params::BendersSeqParam)
        new(oracle, seq_type, params)
    end

    function RootNodePreprocessing(oracle::AbstractOracle; params::BendersSeqParam = BendersSeqParam())
        new(oracle, BendersSeq, params)
    end

    function RootNodePreprocessing(data::Data; params::BendersSeqParam = BendersSeqParam())
        new(ClassicalOracle(data), BendersSeq, params)
    end
end

function root_node_processing!(data::Data, master::AbstractMaster, root_preprocessing::RootNodePreprocessing)
    root_param = deepcopy(root_preprocessing.params)

    undo = relax_integrality(master.model)
    
    root_node_time = @elapsed begin
        BendersRootSeq = root_preprocessing.seq_type(data, master, root_preprocessing.oracle; param=root_param)
        solve!(BendersRootSeq)
    end
    
    undo()
    
    return root_node_time
end

struct EmptyCallbackParam <: AbstractCallbackParam
end

abstract type AbstractLazyCallback end
struct LazyCallback <: AbstractLazyCallback
    params::AbstractCallbackParam
    oracle::AbstractOracle
    
    function LazyCallback(; params=EmptyCallbackParam(), oracle)
        new(params, oracle)
    end
    
    function LazyCallback(data; params=EmptyCallbackParam())
        new(params, ClassicalOracle(data))
    end
end

Base.@kwdef struct UserCallbackParam <: AbstractCallbackParam
    frequency::Int = 50
    node_count::Int = -1
    depth::Int = -1
end

abstract type AbstractUserCallback end
struct UserCallback <: AbstractUserCallback
    params::UserCallbackParam
    oracle::AbstractOracle
    
    function UserCallback(; params=UserCallbackParam(), oracle)
        new(params, oracle)
    end
    
    function UserCallback(data; params=UserCallbackParam())
        new(params, ClassicalOracle(data))
    end
end

function lazy_callback(cb_data, master_model::Model, log::BendersBnBLog, callback::LazyCallback)
    status = JuMP.callback_node_status(cb_data, master_model)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER

        state = BendersBnBState()
        if solver_name(master_model) == "CPLEX"
            n_count = Ref{CPXINT}()
            ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)
            state.node = n_count[]
        end

        state.values[:x] = JuMP.callback_value.(cb_data, master_model[:x])
        state.values[:t] = JuMP.callback_value.(cb_data, master_model[:t])

        state.oracle_time = @elapsed begin
            state.is_in_L, hyperplanes, state.f_x = generate_cuts(callback.oracle, state.values[:x], state.values[:t])
            cuts = !state.is_in_L ? hyperplanes_to_expression(master_model, hyperplanes, master_model[:x], master_model[:t]) : []
        end

        if !isempty(cuts)
            for cut in cuts
                cut_constraint = @build_constraint(0 >= cut)
                MOI.submit(master_model, MOI.LazyConstraint(cb_data), cut_constraint)
                state.num_cuts += 1
            end
        end
        record_node!(log, state, true)
    end
end

function user_callback(cb_data, master_model::Model, log::BendersBnBLog, callback::UserCallback)
    status = JuMP.callback_node_status(cb_data, master_model)
    frequency = callback.params.frequency
    node_count = callback.params.node_count == -1 ? callback.params.node_count : Inf
    depth = callback.params.depth == -1 ? callback.params.depth : 0
    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        log.num_of_fraction_node += 1
        if log.num_of_fraction_node >= frequency 
            if solver_name(master_model) == "CPLEX"
                n_count = Ref{CPXINT}()
                ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)
                depth = Ref{CPXINT}()
                ret2 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODEDEPTH, depth)
                if n_count[] <= node_count && depth[] >= depth
                    log.num_of_fraction_node = 0
                    state = BendersBnBState()
                    state.values[:x] = JuMP.callback_value.(cb_data, master_model[:x])
                    state.values[:t] = JuMP.callback_value.(cb_data, master_model[:t])

                    state.oracle_time = @elapsed begin
                        state.is_in_L, hyperplanes, state.f_x = generate_cuts(callback.oracle, state.values[:x], state.values[:t])

                        cuts = !state.is_in_L ? hyperplanes_to_expression(master_model, hyperplanes, master_model[:x], master_model[:t]) : []
                    end
    
                    if !isempty(cuts)
                        for cut in cuts
                            cut_constraint = @build_constraint(0 >= cut)
                            MOI.submit(master_model, MOI.UserCut(cb_data), cut_constraint)
                            state.num_cuts += 1
                        end
                    end
                    record_node!(log, state, false)
                end
            else
                @warn "node_count and depth are not supported for this solver"
                log.num_of_fraction_node = 0
                state = BendersBnBState()
                state.values[:x] = JuMP.callback_value.(cb_data, master_model[:x])
                state.values[:t] = JuMP.callback_value.(cb_data, master_model[:t])
                state.oracle_time = @elapsed begin
                    state.is_in_L, hyperplanes, state.f_x = generate_cuts(callback.oracle, state.values[:x], state.values[:t])
                    cuts = !state.is_in_L ? hyperplanes_to_expression(master_model, hyperplanes, master_model[:x], master_model[:t]) : []
                end
                if !isempty(cuts)
                    for cut in cuts
                        cut_constraint = @build_constraint(0 >= cut)
                        MOI.submit(master_model, MOI.UserCut(cb_data), cut_constraint)
                        state.num_cuts += 1
                    end
                end
                record_node!(log, state, false)
            end
        end
    end
end