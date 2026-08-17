
"""
    BendersSeq <: AbstractBendersSeq

Sequential Benders decomposition environment.

`BendersSeq` implements the classical Benders cutting-plane algorithm by repeatedly solving the master problem, evaluating the oracle at the current candidate solution, adding the resulting Benders cuts to the master problem, and continuing until a termination criterion is satisfied.

# Fields
- `master::AbstractMaster`: Master problem formulation
- `oracle::AbstractOracle`: Oracle used for subproblem evaluation and cut generation.
- `param::BendersSeqParam`: Parameters controlling the algorithm, such as time limit, gap tolerance, and verbosity
- `obj_value::Float64`: Objective value of the best solution found
- `termination_status::TerminationStatus`: Status of the algorithm upon termination

# Examples
```julia
master = Master(data; model = update_master_model!)
oracle = ClassicalOracle(data, master; model = update_sub_model!)
env = BendersSeq(master, oracle)  # Use default parameters
df = solve!(env)
```

See also: [`BendersSeqInOut`](@ref), [`SpecializedBendersSeq`](@ref), [`BendersBnB`](@ref)
"""
mutable struct BendersSeq <: AbstractBendersSeq
    master::AbstractMaster
    oracle::AbstractOracle

    param::BendersSeqParam

    preprocessing::AbstractBendersPreprocessing

    # result
    obj_value::Float64
    termination_status::TerminationStatus

    function BendersSeq(master::AbstractMaster, oracle::AbstractOracle; param::BendersSeqParam = BendersSeqParam(), preprocessing::AbstractBendersPreprocessing = NoPreprocessing())
        new(master, oracle, param, preprocessing, Inf, NotSolved())
    end
end

"""
    solve!(env::BendersSeq; iter_prefix = "") -> DataFrame

Execute the sequential Benders decomposition algorithm.

At each iteration, the method solves the master problem, queries the oracle
at the resulting candidate solution, updates the bounds and termination state,
and adds any generated Benders cuts to the master problem. Iterations continue
until the problem is solved or a termination condition is reached.

The iteration history is returned as a `DataFrame`.

# Arguments
- `env::BendersSeq`: The configured sequential Benders environment
- `iter_prefix::String`: Optional prefix for iteration logging (default: "")

# Returns
- `DataFrame`: A log of iterations containing lower bounds, upper bounds, gaps, and timing information

# Termination
The algorithm terminates when one of the following occurs:
- the candidate solution is in `L`
- the prescribed time limit is reached;
- the prescribed gap tolerance is satisfied;
- the master problem becomes infeasible or encounters an unexpected status.

The environment's `obj_value` and `termination_status` fields are updated
before returning.
"""
function solve!(env::BendersSeq; iter_prefix = "") 
    log = BendersSeqLog()
    param = env.param
    try    
        # Apply preprocessing
        log.preprocessing_time = preprocess!(env.master, env.preprocessing)

        if param.time_limit <= log.preprocessing_time
            throw(TimeLimitException("Time limit reached before BnB starts, please increase the time limit."))
        end
        
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
                        throw(UnexpectedModelStatusException("BendersSeq: master $(termination_status(env.master.model))"))
                    end
                end
                
                # Execute oracle
                state.oracle_time = @elapsed begin
                    state.is_in_L, hyperplanes, state.f_x = generate_cuts(env.oracle, state.values[:x], state.values[:t]; time_limit = get_sec_remaining(log, param))
                    
                    cuts = !state.is_in_L ? hyperplanes_to_expression(env.master.model, hyperplanes, env.master.x, env.master.t) : []
                
                    if state.f_x !== NaN
                        update_upper_bound_and_gap!(state, log, (f_x, x) -> env.master.c_t' * f_x + env.master.c_x' * x)
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
        @warn e.msg
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
