# ----------------------------------------------------------------------------
# Callback infrastructure
# ----------------------------------------------------------------------------

include("callback/callback.jl") # must be included first
 
"""
    BendersBnB <: AbstractBendersBnB

Branch-and-Bound implementation of Benders decomposition.

`BendersBnB` integrates Benders cut generation into the MIP solver's branch-and-bound process through lazy-constraint and optional user-cut callbacks. Optional root-node preprocessing can be applied before the branch-and-bound search begins.

# Fields
- `master::AbstractMaster`: Master problem.
- `param::BendersBnBParam`: Parameters controlling algorithm, such as the time limit, relative gap tolerance, and verbosity.
- `preprocessing::AbstractBendersPreprocessing`: Configuration for root-node preprocessing.
- `lazy_callback::AbstractLazyCallback`: Configuration for lazy-constraint callback
- `user_callback::AbstractUserCallback`: Configuration for user-cut callback
- `obj_value::Float64`: Objective value of the best solution found
- `termination_status::TerminationStatus`: Status of the algorithm at termination

# Examples
```julia
master = Master(data; model = update_master_model!)
oracle = ClassicalOracle(data, master; model = update_sub_model!)
env = BendersBnB(master, oracle)  # Use default setting with no root node preprocessing and no user callback
result = solve!(env)
```

See also: [`BendersSeq`](@ref)
"""
mutable struct BendersBnB <: AbstractBendersBnB
    master::AbstractMaster 

    param::BendersBnBParam 

    preprocessing::AbstractBendersPreprocessing
    lazy_callback::AbstractLazyCallback
    user_callback::AbstractUserCallback

    obj_value::Float64 
    termination_status::TerminationStatus 

    function BendersBnB(master::AbstractMaster, oracle::AbstractTypicalOracle; param::BendersBnBParam = BendersBnBParam())
        
        preprocessing = NoPreprocessing()
        lazy_callback = LazyCallback(oracle)
        user_callback = NoUserCallback()
        
        new(master, param, preprocessing, lazy_callback, user_callback, Inf, NotSolved())
    end

    function BendersBnB(master::AbstractMaster, preprocessing::AbstractBendersPreprocessing, lazy_callback::AbstractLazyCallback, user_callback::AbstractUserCallback; param::BendersBnBParam = BendersBnBParam())
        
        new(master, param, preprocessing, lazy_callback, user_callback, Inf, NotSolved())
    end
end


"""
    solve!(env::BendersBnB) -> DataFrame

Execute the branch-and-bound Benders decomposition algorithm.

The method applies the configured root-node preprocessing, installs the configured lazy-constraint and user-cut callbacks, sets the solver parameters, and solves the master problem through the MIP solver's branch-and-bound procedure.

The returned `DataFrame` contains summary statistics for the solve.

# Arguments
- `env::BendersBnB`: The configured Benders Branch-and-Bound environment

# Returns
A one-row `DataFrame` containing the solve statistics, including the objective value, objective bound, elapsed time, relative gap, and cut-generation counts. The objective value is Inf when no feasible solution is available.

# Termination
The environment's `obj_value` and `termination_status` fields are updated before returning. A time limit or unexpected solver status is handled and reflected in the returned environment and solve statistics.
"""
function solve!(env::BendersBnB) 
    log = BendersBnBLog()
    param = env.param
    try 
        log.start_time = time()
        
        # Apply root node preprocessing if specified
        log.preprocessing_time = preprocess!(env.master, env.preprocessing)

        # Set up lazy callback
        function lazy_callback_wrapper(cb_data)
            lazy_callback(cb_data, env.master, log, env.param, env.lazy_callback)
        end
        set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback_wrapper)
        
        # Set up user callback if specified
        if !isa(env.user_callback, NoUserCallback)
            function user_callback_wrapper(cb_data)
                user_callback(cb_data, env.master, log, env.param, env.user_callback)
            end
            set_attribute(env.master.model, MOI.UserCutCallback(), user_callback_wrapper)
        end
        
        # Configure solver parameters
        if param.time_limit <= log.preprocessing_time
            throw(TimeLimitException("Time limit reached before BnB starts, please increase the time limit."))
        end
        set_time_limit_sec(env.master.model, param.time_limit - log.preprocessing_time)
        set_optimizer_attribute(env.master.model, MOI.Silent(), !param.verbose)
        set_optimizer_attribute(env.master.model, MOI.RelativeGapTolerance(), param.gap_tolerance)
        
        # Solve the master problem
        JuMP.optimize!(env.master.model)
        
        log.total_time = time() - log.start_time

        # Process termination status
        status = termination_status(env.master.model)
        if status == MOI.OPTIMAL
            env.termination_status = Optimal()
            env.obj_value = JuMP.objective_value(env.master.model)
        elseif status == MOI.TIME_LIMIT
            env.termination_status = TimeLimit()
            env.obj_value = has_values(env.master.model) ? JuMP.objective_value(env.master.model) : Inf
        else
            throw(UnexpectedModelStatusException("BendersBnB: master $(status)"))
        end
       
        df = to_dataframe(env, log)

        if param.verbose 
            @info df
        end
        
        return df
    catch e
        @warn e.msg
        if typeof(e) <: TimeLimitException
            env.termination_status = TimeLimit()
            env.obj_value = has_values(env.master.model) ? JuMP.objective_value(env.master.model) : Inf
        elseif typeof(e) <: UnexpectedModelStatusException
            env.termination_status = InfeasibleOrNumericalIssue()
        else
            rethrow()  
        end
        if env.param.verbose
            println("Terminated with $(env.termination_status)")
        end
        return to_dataframe(env, log)
    end
end
