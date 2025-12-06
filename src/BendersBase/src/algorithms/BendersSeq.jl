export BendersSeq, solve!

"""
    BendersSeq <: AbstractBendersSeq

Sequential Benders decomposition algorithm using a cutting-plane method.

This is the basic Benders decomposition implementation that iteratively solves the master problem, generates Benders cuts from the oracle, and refines the master problem until convergence or a termination criterion is met.

# Fields
- `data::Data`: Problem data containing dimensions, cost vectors, and problem-specific information
- `master::AbstractMaster`: Master problem formulation
- `oracle::AbstractOracle`: Oracle for subproblem solving and cut generation
- `param::BendersSeqParam`: Parameters controlling algorithm behavior (time limit, gap tolerance, verbosity, etc.)
- `obj_value::Float64`: Objective value of the best solution found
- `termination_status::TerminationStatus`: Status of the algorithm upon termination

# Constructors
```julia
BendersSeq(data, master::AbstractMaster, oracle::AbstractOracle; param::BendersSeqParam = BendersSeqParam())
BendersSeq(data; param::BendersSeqParam = BendersSeqParam())  # Uses default Master and ClassicalOracle
```

# Examples
```julia
data = Data(...)
algorithm = BendersSeq(data)  # Use default parameters
df = solve!(algorithm)
```

See also: [`BendersSeqInOut`](@ref), [`SpecializedBendersSeq`](@ref), [`BendersBnB`](@ref)
"""
mutable struct BendersSeq <: AbstractBendersSeq
    data::Data
    master::AbstractMaster
    oracle::AbstractOracle

    param::BendersSeqParam

    # result
    obj_value::Float64
    termination_status::TerminationStatus

    function BendersSeq(data, master::AbstractMaster, oracle::AbstractOracle; param::BendersSeqParam = BendersSeqParam())
        new(data, master, oracle, param, Inf, NotSolved())
    end

    function BendersSeq(data; param::BendersSeqParam = BendersSeqParam())
        new(data, Master(data), ClassicalOracle(data), param, Inf, NotSolved())
    end

    function BendersSeq(master::AbstractMaster, oracle::AbstractOracle; param::BendersSeqParam = BendersSeqParam())

        # Initialize data object
        dim_x = length(master.x)
        obj = objective_function(master.model)
        c_x = [coefficient(obj, master.x[i]) for i in 1:dim_x]
        
        dim_t = length(master.t)
        c_t = [coefficient(obj, master.t[i]) for i in 1:dim_t]
        
        data = Data(dim_x, dim_t, EmptyData(), c_x, c_t)

        new(data, master, oracle, param, Inf, NotSolved())
    end
end

"""
    solve!(env::BendersSeq; iter_prefix = "") -> DataFrame

Execute the sequential Benders decomposition algorithm.

This function implements the core Benders cutting-plane method: repeatedly solving the master problem, querying the oracle for Benders cuts, and adding cuts to refine the master problem.

# Arguments
- `env::BendersSeq`: The configured Benders algorithm environment
- `iter_prefix::String`: Optional prefix for iteration logging (default: "")

# Returns
- `DataFrame`: A log of iterations containing lower bounds, upper bounds, gaps, and timing information

# Algorithm Steps
1. Solve the master problem to obtain candidate solution (x, t)
2. Query the oracle to check feasibility and generate Benders cuts
3. Update bounds and check termination criteria
4. Add generated cuts to the master problem
5. Repeat until convergence or termination

# Termination Criteria
- Optimal solution found (point is in L)
- Time limit reached
- Gap tolerance met
- Master problem becomes infeasible

# Examples
```julia
algorithm = BendersSeq(data)
log_df = solve!(algorithm)
println("Objective: ", algorithm.obj_value)
```
"""
function solve!(env::BendersSeq; iter_prefix = "") 
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
                        update_upper_bound_and_gap!(state, log, (f_x, x) -> env.data.c_t' * f_x + env.data.c_x' * x)
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
end

