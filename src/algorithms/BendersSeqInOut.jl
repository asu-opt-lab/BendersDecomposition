export BendersSeqInOut
mutable struct BendersSeqInOut <: AbstractBendersSeq
    data::Data
    master::AbstractMaster
    oracle::AbstractOracle

    param::BendersSeqInOutParam 

    # result
    obj_value::Float64
    termination_status::TerminationStatus

    function BendersSeqInOut(data, master::AbstractMaster, oracle::AbstractOracle; param::BendersSeqInOutParam = BendersSeqInOutParam()) 
        new(data, master, oracle, param, Inf, NotSolved())
    end

    function BendersSeqInOut(data; param::BendersSeqInOutParam = BendersSeqInOutParam())
        new(data, Master(data), ClassicalOracle(data), param, Inf, NotSolved())
    end
end
"""
Run BendersSeqInOut
"""
# param should have stabilizing_x
function solve!(env::BendersSeqInOut)
    param = env.param
    log = BendersSeqLog()
    try    
        state = BendersSeqState()
        stabilizing_x = param.stabilizing_x
        α = param.α
        λ = param.λ
        # relax_integrality(env.master.model)
        kelley_mode = false
        
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
                        throw(ErrorException("master termination status: $(termination_status(env.master.model))"))
                        # if infeasible, then the milp is infeasible
                    end
                end
                
                # perturb point
                stabilizing_x = α * stabilizing_x + (1 - α) * state.values[:x]
                intermediate_x = λ * state.values[:x] + (1 - λ) * stabilizing_x

                # Execute oracle
                state.oracle_time = @elapsed begin
                    state.is_in_L, hyperplanes, state.f_x = generate_cuts(env.oracle, intermediate_x, state.values[:t]; time_limit = get_sec_remaining(log, param)) # only modify the knapsack cut and classical cut

                    if kelley_mode 
                        if state.f_x != NaN
                            update_upper_bound_and_gap!(state, log, (f_x, x) -> env.data.c_t' * f_x + env.data.c_x' * x)
                        end
                    else
                        state.is_in_L = false
                    end

                    cuts = !state.is_in_L ? hyperplanes_to_expression(env.master.model, hyperplanes, env.master.model[:x], env.master.model[:t]) : []
                end
            
                # Update state and record information
                record_iteration!(log, state)

                param.verbose && print_iteration_info(state, log)
            end

            # Check termination criteria
            is_terminated(state, log, param) && break
            # @info "state.is_in_L: $(state.is_in_L)"
            # @info "cuts: $cuts"
            # add generated cuts to master
            @constraint(env.master.model, 0 .>= cuts)
            
            # whether to switch kelley mode
            if !kelley_mode && log.n_iter != 0
                check_lb_improvement!(state, log; zero_tol = 1e-8, tol_imprv = 0.05)

                if log.consecutive_no_improvement >= 5
                    # Reset λ to 1 (switch to Kelley's cutting plane)
                    λ = 1.0
                    kelley_mode = true
                    param.verbose && println("Switching to Kelley's cutting plane method (λ = 1.0)")
                end
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
    # even if it terminates in the middle due to time limit, should be able to access the latest x_value via env.iterations[end].x_value
end