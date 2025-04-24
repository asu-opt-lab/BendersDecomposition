export SpecializedBendersSeq

mutable struct SpecializedBendersSeq <: AbstractBendersSeq
    data::Data
    master::AbstractMaster
    oracle::DisjunctiveOracle

    param::SpecializedBendersSeqParam # initially default and add an interface function?

    # result
    obj_value::Float64
    termination_status::TerminationStatus

    function SpecializedBendersSeq(data, master::AbstractMaster, oracle::DisjunctiveOracle; param::SpecializedBendersSeqParam = SpecializedBendersSeqParam()) 
        relax_integrality(master.model)

        oracle.oracle_param.split_index_selection_rule != LargestFractional() && throw(AlgorithmException("SpeicalizedBendersSeq does not admit $(oracle.oracle_param.split_index_selection_rule). Use LargestFractional() instead."))
        oracle.oracle_param.disjunctive_cut_append_rule != DisjunctiveCutsSmallerIndices() && throw(AlgorithmException("SpeicalizedBendersSeq does not admit $(oracle.oracle_param.disjunctive_cut_append_rule). Use DisjunctiveCutsSmallerIndices() instead."))         

        # case where master and oracle has their own attributes and default loop_param and solver_param
        new(data, master, oracle, param, Inf, NotSolved())
    end
end

"""
Run BendersSeq
"""
function solve!(env::SpecializedBendersSeq) 
    log = BendersSeqLog()
    L_param = BendersSeqParam(; time_limit = env.param.time_limit, gap_tolerance = env.param.gap_tolerance, verbose = env.param.verbose)
    L_env = BendersSeq(env.data, env.master, env.oracle.typical_oracles[1]; param = L_param)

    try
        while true
            state = BendersSeqState()
            state.total_time = @elapsed begin
                # Solve linear relaxation
                state.master_time = @elapsed begin
                    solve!(L_env)
                    state.LB, state.values[:x], state.values[:t] = JuMP.objective_value(env.master.model), value.(env.master.model[:x]), value.(env.master.model[:t])
                end
                @debug value.(env.master.model[:x])

                # Check termination criteria
                is_terminated(state, log, env.param) && (record_iteration!(log, state); break)

                # Execute oracle
                state.oracle_time = @elapsed begin
                    state.is_in_L, hyperplanes, state.f_x = generate_cuts(env.oracle, state.values[:x], state.values[:t]; time_limit = get_sec_remaining(log, env.param))
                    cuts = !state.is_in_L ? hyperplanes_to_expression(env.master.model, hyperplanes, env.master.model[:x], env.master.model[:t]) : []
                end

                state.is_in_L == true && throw(InfeasibleOrNumericalIssue("τ=0 at fractional point, possibily numerical issue"))
                
                record_iteration!(log, state)
                env.param.verbose && print_iteration_info(state, log)

                # Add generated cuts to master
                @constraint(env.master.model, 0 .>= cuts)
            end
        end
        env.termination_status = Optimal()
        env.obj_value = log.iterations[end].LB
        
        return to_dataframe(log)
    catch e
        if typeof(e) <: TimeLimitException
            env.termination_status = TimeLimit()
            env.obj_value = log.iterations[end].LB
        elseif typeof(e) <: UnexpectedModelStatusException
            env.termination_status = InfeasibleOrNumericalIssue()
        else
            rethrow()  
        end
        return to_dataframe(log)
    end
# even if it terminates in the middle due to time limit, should be able to access the latest x_value via env.iterations[end].values[:x]
end
