export BendersSeqInOut

"""
    BendersSeqInOut <: AbstractBendersSeq

Sequential Benders decomposition with In-Out stabilization technique.

This variant of Benders decomposition uses stabilization to improve convergence by maintaining a stabilizing point and perturbing query points. The algorithm can switch to Kelley's cutting-plane method if progress stalls.

# Fields
- `master::AbstractMaster`: Master module
- `oracle::AbstractOracle`: Oracle module for cut generation
- `param::BendersSeqInOutParam`: Parameters controlling algorithm behavior (includes stabilization parameters α, λ, and stabilizing_x)
- `obj_value::Float64`: Objective value of the best solution found
- `termination_status::TerminationStatus`: Status of the algorithm upon termination

# Constructors
```julia
BendersSeqInOut(master::AbstractMaster, oracle::AbstractOracle; param::BendersSeqInOutParam = BendersSeqInOutParam())
```

# Stabilization Parameters
The stabilization technique requires three parameters (specified in `BendersSeqInOutParam`):
- `α`: Weight for updating the stabilizing point
- `λ`: Weight for perturbing the query point
- `stabilizing_x`: Initial stabilizing point

# Examples
```julia
master = Master(data; customize = customize_master_model!)
oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
param = BendersSeqInOutParam(α = 0.8, λ = 0.5, stabilizing_x = zeros(master.dim_x))
env = BendersSeqInOut(master, oracle; param = param)
df = solve!(env)
```

See also: [`BendersSeq`](@ref), [`SpecializedBendersSeq`](@ref)
"""
mutable struct BendersSeqInOut <: AbstractBendersSeq
    master::AbstractMaster
    oracle::AbstractOracle

    param::BendersSeqInOutParam 

    # result
    obj_value::Float64
    termination_status::TerminationStatus

    function BendersSeqInOut(master::AbstractMaster, oracle::AbstractOracle; param::BendersSeqInOutParam = BendersSeqInOutParam())

        new(master, oracle, param, Inf, NotSolved())
    end
end
"""
    solve!(env::BendersSeqInOut) -> DataFrame

Execute the sequential Benders decomposition with In-Out stabilization.

This function implements a stabilized Benders cutting-plane method that uses a perturbed query point
to improve convergence. If the lower bound stagnates for multiple iterations, the algorithm switches
to Kelley's cutting-plane method (setting λ = 1.0).

# Arguments
- `env::BendersSeqInOut`: The configured Benders In-Out algorithm environment

# Returns
- `DataFrame`: A log of iterations containing lower bounds, upper bounds, gaps, and timing information

# Algorithm Steps
1. Solve the master problem to obtain candidate solution (x, t)
2. Update the stabilizing point: `stabilizing_x = α * stabilizing_x + (1 - α) * x`
3. Generate perturbed query point: `intermediate_x = λ * x + (1 - λ) * stabilizing_x`
4. Query the oracle at the perturbed point to generate Benders cuts
5. Update bounds and check termination criteria
6. Add generated cuts to the master problem
7. Check if switching to Kelley mode is needed (if no improvement for 5 consecutive iterations)
8. Repeat until convergence or termination

# Termination Criteria
- Optimal solution found (point is in L after switching to Kelley mode)
- Time limit reached
- Gap tolerance met
- Master problem becomes infeasible
"""
function solve!(env::BendersSeqInOut)
    param = env.param
    log = BendersSeqLog()
    try    
        state = BendersSeqState()
        stabilizing_x = param.stabilizing_x
        α = param.α
        λ = param.λ
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
                        state.values[:x] = JuMP.value.(env.master.x)
                        state.values[:t] = JuMP.value.(env.master.t)
                    elseif termination_status(env.master.model) == TIME_LIMIT
                        throw(TimeLimitException("Time limit reached during master solving"))
                    else
                        throw(ErrorException("master termination status: $(termination_status(env.master.model))"))
                    end
                end
                
                # perturb point
                stabilizing_x = α * stabilizing_x + (1 - α) * state.values[:x]
                intermediate_x = λ * state.values[:x] + (1 - λ) * stabilizing_x

                # Execute oracle
                state.oracle_time = @elapsed begin
                    state.is_in_L, hyperplanes, state.f_x = generate_cuts(env.oracle, intermediate_x, state.values[:t]; time_limit = get_sec_remaining(log, param))

                    if kelley_mode 
                        if state.f_x != NaN
                            update_upper_bound_and_gap!(state, log, (f_x, x) -> env.master.c_t' * f_x + env.master.c_x' * x)
                        end
                    else
                        state.is_in_L = false
                    end

                    cuts = !state.is_in_L ? hyperplanes_to_expression(env.master.model, hyperplanes, env.master.x, env.master.t) : []
                end
            
                # Update state and record information
                record_iteration!(log, state)
            end

            param.verbose && print_iteration_info(state, log)

            # Check termination criteria
            is_terminated(state, log, param) && break

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
            @warn e.msg
        else
            rethrow()  
        end
        if env.param.verbose
            println("Terminated with $(env.termination_status)")
        end
        return to_dataframe(log)
    end
end