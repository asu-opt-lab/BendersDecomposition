export SpecializedBendersSeq, solve!

include(joinpath(dirname(dirname(@__DIR__)), "example", "uflp", "model.jl"))

mutable struct SpecializedBendersSeq <: AbstractBendersSeq
    data::Data
    master::AbstractMaster
    typical_oracle::AbstractOracle
    oracle::AbstractOracle

    param::SpecializedBendersSeqParam # initially default and add an interface function?

    # result
    obj_value::Float64
    termination_status::TerminationStatus

    function SpecializedBendersSeq(data, master::AbstractMaster, typical_oracle::AbstractOracle, oracle::DisjunctiveOracle; param::SpecializedBendersSeqParam = SpecializedBendersSeqParam()) 
        relax_integrality(master.model)

        typeof(typical_oracle) == DisjunctiveOracle && throw(UndefError("assign TypicalOracle to typical_oracle instead of $(typeof(typical_oracle))"))
        oracle.oracle_param.split_index_selection_rule != LargestFractional() && throw(UndefError("assign LargestFractional() to split_index_selection_rule instead of $(oracle.oracle_param.split_index_selection_rule)"))
        oracle.oracle_param.disjunctive_cut_append_rule != DisjunctiveCutsSmallerIndices() && throw(UndefError("assign DisjunctiveCutsSmallerIndices() to disjunctive_cut_append_rule instead of $(oracle.oracle_param.disjunctive_cut_append_rule)"))         

        # case where master and oracle has their own attributes and default loop_param and solver_param
        new(data, master, typical_oracle, oracle, param, Inf, NotSolved())
    end

    function SpecializedBendersSeq(data; param::SpecializedBendersSeqParam = SpecializedBendersSeqParam())
        master = Master(data)
        update_model!(master, data)
        relax_integrality(master.model)

        typical_oracles = [ClassicalOracle(data); ClassicalOracle(data)] 
        map(k -> update_model!(typical_oracles[k], data), 1:2)

        disjunctive_oracle = DisjunctiveOracle(data, typical_oracles)

        set_parameter!(disjunctive_oracle, "split_index_selection_rule", LargestFractional())
        set_parameter!(disjunctive_oracle, "disjunctive_cut_append_rule", DisjunctiveCutsSmallerIndices())

        typical_oracle = ClassicalOracle(data)
        update_model!(typical_oracle, data)
        
        # case where master and oracle has their own attributes and default loop_param and solver_param
        new(data, master, typical_oracle, disjunctive_oracle, param, Inf, NotSolved())
    end

    function SpecializedBendersSeq(data, master::AbstractMaster, typical_oracle::AbstractOracle, oracle::AbstractOracle; param::SpecializedBendersSeqParam = SpecializedBendersSeqParam())
        throw(UndefError("assign DisjunctiveOracle to oracle instead of $(typeof(oracle))"))
    end
end

"""
Run BendersSeq
"""
function solve!(env::SpecializedBendersSeq) 
    log = BendersSeqLog()
    L_param = BendersSeqParam(; time_limit = env.param.time_limit, gap_tolerance = env.param.gap_tolerance, verbose = env.param.verbose)
    L_env = BendersSeq(env.data, env.master, env.typical_oracle; param = L_param)

    try
        while true
            state = BendersSeqState()
            state.total_time = @elapsed begin
                # Solve linear relaxation
                state.master_time = @elapsed begin
                    solve!(L_env)
                    state.LB, state.values[:x], state.values[:t] = JuMP.objective_value(env.master.model), value.(env.master.model[:x]), value.(env.master.model[:t])
                end
                println(value.(env.master.model[:x]))
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
