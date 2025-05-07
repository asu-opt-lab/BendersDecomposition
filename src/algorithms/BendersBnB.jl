export BendersBnB

include("callback/preprocessing.jl") # must be included first
include("callback/callback.jl") # must be included first
 
"""
    BendersBnB <: AbstractBendersCallback

Branch-and-Bound implementation of Benders decomposition algorithm.

This implementation uses callbacks to efficiently generate Benders cuts during the branch-and-bound process,
avoiding the need to repeatedly solve the entire master problem.

# Fields
- `data::Data`: Problem data containing dimensions, cost vectors, and problem-specific information
- `master::AbstractMaster`: Master problem formulation
- `param::BendersBnBParam`: Parameters controlling algorithm behavior
- `root_preprocessing::AbstractRootNodePreprocessing`: Configuration for preprocessing at the root node
- `lazy_callback::AbstractLazyCallback`: Configuration for lazy constraint callbacks
- `user_callback::AbstractUserCallback`: Configuration for user cut callbacks
- `obj_value::Float64`: Objective value of the best solution found
- `termination_status::TerminationStatus`: Status of the algorithm upon termination

# Examples
```julia
data = Data(...)  # Create problem data
algorithm = BendersBnB(data)  # Use default parameters
obj_value, solve_time = solve!(algorithm)
```
"""
mutable struct BendersBnB <: AbstractBendersCallback
    data::Data
    master::AbstractMaster 

    param::BendersBnBParam 

    root_preprocessing::AbstractRootNodePreprocessing
    lazy_callback::AbstractLazyCallback
    user_callback::AbstractUserCallback

    obj_value::Float64 
    termination_status::TerminationStatus 

    function BendersBnB(data; param::BendersBnBParam = BendersBnBParam())
        new(data, Master(data), param, RootNodePreprocessing(data), LazyCallback(data), UserCallback(data), Inf, NotSolved())
    end

    function BendersBnB(data, master::AbstractMaster, root_preprocessing::AbstractRootNodePreprocessing, lazy_callback::AbstractLazyCallback, user_callback::AbstractUserCallback; param::BendersBnBParam = BendersBnBParam())
        new(data, master, param, root_preprocessing, lazy_callback, user_callback, Inf, NotSolved())
    end
end

"""
    solve!(env::BendersBnB) -> Tuple{Float64, Float64}

Execute the branch-and-bound Benders decomposition algorithm.

This function configures callbacks, solves the master problem with the callback-based 
cutting plane approach, and processes the results.

# Arguments
- `env::BendersBnB`: The configured Benders Branch-and-Bound algorithm environment

# Returns
- `Tuple{Float64, Float64}`: A tuple containing (objective_value, elapsed_time)
  - The objective value will be Inf if no feasible solution is found
  - Elapsed time is measured in seconds

# Algorithm Steps
1. Apply root node preprocessing if specified (relaxing integrality constraints)
2. Configure lazy and user callbacks for the branch-and-bound process
3. Set solver parameters (time limit, verbosity, and gap tolerance)
4. Solve the master problem with callbacks
5. Process termination status and objective value
6. Return results and execution statistics
"""
function solve!(env::BendersBnB) 
    log = BendersBnBLog()
    param = env.param
    log.start_time = time()
    
    # Apply root node preprocessing if specified
    root_node_time = 0.0
    if isa(env.root_preprocessing, RootNodePreprocessing)
        root_node_time = root_node_processing!(env.data, env.master, env.root_preprocessing)
    end
    
    # Apply disjunctive root node preprocessing if specified
    if param.disjunctive_root_process
        # Update root_prepreocessing params
        env.root_preprocessing.params.time_limit -= root_node_time
        env.root_preprocessing.oracle = env.user_callback.oracle

        disjunctive_root_node_time = root_node_processing!(env.data, env.master, env.root_preprocessing)
        root_node_time += disjunctive_root_node_time
    end
    
    # Set up lazy callback
    function lazy_callback_wrapper(cb_data)
        lazy_callback(cb_data, env.master.model, log, env.param, env.lazy_callback)
    end
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback_wrapper)
    
    # Set up user callback if specified
    if !isa(env.user_callback, NoUserCallback)
        function user_callback_wrapper(cb_data)
            user_callback(cb_data, env.master.model, log, env.param, env.user_callback)
        end
        set_attribute(env.master.model, MOI.UserCutCallback(), user_callback_wrapper)
    end
    
    # Configure solver parameters
    if param.time_limit <= root_node_time
        throw(TimeLimitException("Time limit reached before BnB starts, please increase the time limit."))
    end
    set_time_limit_sec(env.master.model, param.time_limit - root_node_time)
    set_optimizer_attribute(env.master.model, MOI.Silent(), !param.verbose)
    set_optimizer_attribute(env.master.model, MOI.RelativeGapTolerance(), param.gap_tolerance)
    
    # Solve the master problem
    JuMP.optimize!(env.master.model)
    
    # Process termination status
    status = termination_status(env.master.model)
    if status == MOI.OPTIMAL
        env.termination_status = Optimal()
        env.obj_value = JuMP.objective_value(env.master.model)
    elseif status == MOI.TIME_LIMIT
        env.termination_status = TimeLimit()
        env.obj_value = has_values(env.master.model) ? JuMP.objective_value(env.master.model) : Inf
    else
        env.termination_status = InfeasibleOrNumericalIssue()
        env.obj_value = Inf
    end
    
    elapsed_time = time() - log.start_time
    
    # Print summary if verbose mode is enabled
    if param.verbose 
        @info "Node count: $(JuMP.node_count(env.master.model))"
        @info "Root processing time: $(root_node_time) "
        @info "Elapsed time: $(elapsed_time)"
        @info "Objective bound: $(JuMP.objective_bound(env.master.model))"
        @info "Objective value: $(env.obj_value)"
        @info "Relative gap: $(JuMP.relative_gap(env.master.model))"
        @info "Lazy cuts added: $(log.n_lazy_cuts)"
        if typeof(env.user_callback.oracle) <: AbstractDisjunctiveOracle
            @info "Disjunctive cuts added: $(length(env.user_callback.oracle.disjunctiveCuts))"
            env.user_callback.oracle.oracle_param.add_benders_cuts_to_master != 0 && @info "Byproduct Benders cuts added: $(log.n_user_cuts - length(env.user_callback.oracle.disjunctiveCuts))"
        end
    end
    
    return env.obj_value, elapsed_time
end



