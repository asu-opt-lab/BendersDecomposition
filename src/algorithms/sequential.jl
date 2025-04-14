"""
Iteration info for Sequential Benders run
"""
mutable struct BendersSeqState <: AbstractBendersSeqState
    master_time::Float64
    oracle_time::Float64
    total_time::Float64
    x_value::Vector{Float64} # should we store this info? move it to log? # or we may keep just last iteration of state. # param for recording states.
    t_value::Vector{Float64} # move it to log?
    f_x::Vector{Float64}
    is_in_L::Bool
    LB::Float64
    UB::Float64
    gap::Float64
   
    # Constructor with specified values
    function BendersSeqState()
        new(0.0, 0.0, 0.0, Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), false, -Inf, Inf, 100.0)
    end
end

"""
Log for Sequential Benders run
"""
mutable struct BendersSeqLog <: AbstractBendersSeqLog
    n_iter::Int
    iterations::Vector{BendersSeqState}
    start_time::Float64
    consecutive_no_improvement::Int
    
    function BendersSeqLog()
        new(0, Vector{BendersSeqState}(), time(), 0)
    end
end

"""
Run BendersSeq
"""
function solve!(env::BendersEnv, ::Seq, params::BendersParams) 
    log = BendersSeqLog()
    try    
        while true
            state = BendersSeqState()
            state.total_time = @elapsed begin
                # Solve master problem
                state.master_time = @elapsed begin
                    set_time_limit_sec(env.master.model, get_sec_remaining(log, params))
                    optimize!(env.master.model)
                    if is_solved_and_feasible(env.master.model; allow_local = false, dual = false)
                        state.LB = JuMP.objective_value(env.master.model)
                        state.x_value = value.(env.master.model[:x])
                        state.t_value = value.(env.master.model[:t])
                    elseif termination_status(env.master.model) == TIME_LIMIT
                        throw(TimeLimitException("Time limit reached during master solving"))
                    else 
                        throw(UnexpectedModelStatusException("Master: $(termination_status(env.master.model))"))
                        # if infeasible, then the milp is infeasible
                    end
                end
                
                # Execute oracle
                state.oracle_time = @elapsed begin
                    state.is_in_L, hyperplanes, state.f_x = generate_cuts(env.oracle, state.x_value, state.t_value; time_limit = get_sec_remaining(log, params))
                    
                    cuts = !state.is_in_L ? hyperplanes_to_expression(env.master.model, hyperplanes, env.master.model[:x], env.master.model[:t]) : []
                
                    if state.f_x !== NaN
                        update_upper_bound_and_gap!(state, (f_x, x) -> env.data.c_t' * f_x + env.data.c_x' * x)
                    end
                end

                # Add generated cuts to master
                @constraint(env.master.model, 0 .>= cuts)

                params.verbose && print_iteration_info(state, log)
            end

            # Update state and record information
            # if params.record_iterations
                record_iteration!(log, state)
            # else
            #     log.iterations[end] = state
            #     log.n_iter += 1
            # end

            # Check termination criteria
            is_terminated(state, log, params) && break
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
    end
# even if it terminates in the middle due to time limit, should be able to access the latest x_value via env.iterations[end].x_value
end

function print_iteration_info(state::BendersSeqState, log::BendersSeqLog)
    @printf("Iter: %4d | LB: %12.4f | UB: %11.4f | Gap: %8.3f%% | Time: (M: %6.2f, S: %6.2f, Total: %6.2f) \n",
           log.n_iter, state.LB, state.UB, state.gap, 
           state.master_time, state.oracle_time, state.total_time)
end

"""
Check termination criteria
"""
function is_terminated(state::BendersSeqState, log::BendersSeqLog, params::BendersParams)
    return state.is_in_L || state.gap <= params.gap_tolerance || get_sec_remaining(log, params) <= 0.0
end
