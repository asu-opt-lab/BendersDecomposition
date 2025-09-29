export DualDecompositionBendersSeq, solve!

include("$(dirname(dirname(@__DIR__)))/example/cflp/oracle.jl")

mutable struct DualDecompositionBendersSeq <: AbstractBendersSeq
    data::Data
    master::AbstractMaster
    typical_oracle::AbstractOracle
    DD_oracle::AbstractPrimalDualOracle

    param::BendersSeqParam # initially default and add an interface function?

    # result
    obj_value::Float64
    termination_status::TerminationStatus

    function DualDecompositionBendersSeq(data, master::AbstractMaster, typical_oracle::AbstractOracle, DD_oracle::AbstractPrimalDualOracle; param::BendersSeqParam = BendersSeqParam()) 
        # case where master and oracle has their own attributes and default loop_param and solver_param
        new(data, master, typical_oracle, DD_oracle, param, Inf, NotSolved())
    end

    function DualDecompositionBendersSeq(data; param::BendersSeqParam = BendersSeqParam())
        # case where master and oracle has their own attributes and default loop_param and solver_param
        new(data, Master(data), ClassicalOracle(data), DualDecomposition(data), param, Inf, NotSolved())
    end
end

"""
Run DualDecompositionBendersSeq
"""
function solve!(env::DualDecompositionBendersSeq; iter_prefix = "") 
    log = BendersSeqLog()
    param = env.param
    try
        while true
            state = BendersSeqState()
            state.total_time = @elapsed begin
                # Solve master problem
                state.master_time = @elapsed begin
                    set_time_limit_sec(env.master.model, get_sec_remaining(log, param))
                    optimize!(env.master.model)
                    if is_solved_and_feasible(env.master.model; allow_local = false, dual = false)
                        state.LB = JuMP.objective_value(env.master.model)
                        state.values[:x] = value.(env.master.model[:x])
                        state.values[:t] = value.(env.master.model[:t])
                    elseif termination_status(env.master.model) == TIME_LIMIT
                        throw(TimeLimitException("Time limit reached during master solving"))
                    else 
                        throw(UnexpectedModelStatusException("BendersSeq: master $(termination_status(env.master.model))"))
                        # if infeasible, then the milp is infeasible
                    end
                end
                
                # Execute dual decomposition
                state.oracle_time = @elapsed begin
                    if length(log.iterations) > env.DD_oracle.oracle_param.LB_stag_consecutive_iter && isapprox(log.iterations[end].LB, log.iterations[end-3].LB)
                        @info "Cplex"
                        state.cplex_iter = log.n_iter + 1
                        state.oracle_time = @elapsed begin
                            state.is_in_L, hyperplanes, state.f_x = generate_cuts(env.typical_oracle, state.values[:x], state.values[:t]; time_limit = get_sec_remaining(log, param))
                            cuts = !state.is_in_L ? hyperplanes_to_expression(env.master.model, hyperplanes, env.master.model[:x], env.master.model[:t]) : []
                        
                            if !isnan(state.f_x[1])
                                update_upper_bound_and_gap!(state, log, (f_x, x) -> env.data.c_t' * f_x + env.data.c_x' * x)
                            end
                        end
                    else
                        @info "Dual Decomposition"
                        t = @elapsed env.DD_oracle.oracle_param.obj_limit = vogel(env.data, state.values[:x]); @info (env.DD_oracle.oracle_param.obj_limit, t)
                        hyperplanes = generate_cuts(env.data, env.DD_oracle, state.values[:x])
                        cuts = hyperplanes_to_expression(env.master.model, hyperplanes, env.master.model[:x], env.master.model[:t])
                        env.DD_oracle.oracle_log.flag_exceed && (state.exceed_iter = log.n_iter + 1; env.DD_oracle.oracle_log.flag_exceed = false)
                    end
                end
                record_iteration!(log, state)
            end
            param.verbose && print_iteration_info(state, log; prefix=iter_prefix)

            # Check termination criteria
            is_terminated(state, log, param) && break

            # Add generated cuts to master
            @constraint(env.master.model, 0.0 .>= cuts)
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
            @warn e.msg
        else
            rethrow()  
        end
        if env.param.verbose
            println("Terminated with $(env.termination_status)")
        end
        return to_dataframe(log)
    end
# even if it terminates in the middle due to time limit, should be able to access the latest x_value via env.iterations[end].values[:x]
end

