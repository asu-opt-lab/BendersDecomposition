export SpecializedBendersSeq, solve!

mutable struct SpecializedBendersSeq <: AbstractBendersSeq
    data::Data
    master::AbstractMaster
    typical_oracle::AbstractOracle
    oracle::AbstractOracle

    param::SpecializedBendersSeqParam # initially default and add an interface function?

    # result
    obj_value::Float64
    termination_status::TerminationStatus

    # Need to handle the following cases when initializing SpecializedBendersSeq:
    # 1. specify typical & specify disjunctive
    # 2. not specify typical & specify disjunctive
    # 3. specify typical & not specify disjunctive
    # 4. not specify typical & not specify disjunctive

    function SpecializedBendersSeq(data, master::AbstractMaster, typical_oracle::AbstractOracle, oracle::DisjunctiveOracle; param::SpecializedBendersSeqParam = SpeializedBendersSeqParam()) 
        relax_integrality(master.model)
        # case where master and oracle has their own attributes and default loop_param and solver_param
        new(data, master, typical_oracle, oracle, param, Inf, NotSolved())
    end

    # function SpecializedBendersSeq(data; param::BendersSeqParam = BendersSeqParam())
    #     relax_integrality(master.model)
    #     # case where master and oracle has their own attributes and default loop_param and solver_param
    #     new(data, Master(data), ClassicalOracle(data), param, Inf, NotSolved())
    # end

    function SpecializedBendersSeq(data, master::AbstractMaster, oracle::AbstractOracle; param::SpecializedBendersSeqParam = SpeializedBendersSeqParam())
        throw(UndefError("assign DisjunctiveOrcale to SpecializedBendersSeq instead of $(typeof(oracle))"))
    end

    function SpecializedBendersSeq(data, master::AbstractMaster, typical_oracle::DisjunctiveOracle, oracle::AbstractOracle; param::SpecializedBendersSeqParam = SpeializedBendersSeqParam())
        throw(UndefError("assign typical oracle instead of $(typeof(oracle))"))
    end
end

"""
Run BendersSeq
"""
function solve!(env::SpecializedBendersSeq) 
    # log = SpecializedBendersSeqLog()
    log = BendersSeqLog()
    L_param = BendersSeqParam(; time_limit = env.param.time_limit, gap_tolerance = env.param.gap_tolerance, verbose = env.param.verbose)
    L_env = BendersSeq(env.data, env.master, env.typical_oracle; param = L_param)

    try
        while true
            # state = SpecializedBendersSeqState()
            state = BendersSeqState()
            state.total_time = @elapsed begin
                # Solve linear relaxation
                state.master_time = @elapsed begin
                    solve!(L_env)
                    state.LB, state.values[:x], state.values[:t] = JuMP.objective_value(env.master.model), value.(env.master.model[:x]), value.(env.master.model[:t])
                end
                
                # Check termination criteria
                is_terminated(state, log, env.param) && (record_iteration!(log, state); break)

                # Execute oracle
                state.oracle_time = @elapsed begin
                    state.is_in_L, hyperplanes, state.f_x = generate_cuts(env.oracle, state.values[:x], state.values[:t]; time_limit = get_sec_remaining(log, env.param))
                    println("is the point in L? $(state.is_in_L)")
                    cuts = !state.is_in_L ? hyperplanes_to_expression(env.master.model, hyperplanes, env.master.model[:x], env.master.model[:t]) : []
                end
                
                record_iteration!(log, state)
                env.param.verbose && print_iteration_info(state, log)
                println(cuts)
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
