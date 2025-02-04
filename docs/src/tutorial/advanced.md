```@meta
EditURL = "advanced.jl"
```

# Advanced Tutorial

For the researchers and industry practitioners, they may need to customize the Benders Decomposition algorithm to solve their own problems.
We provide the `BendersOracle` interface to allow users to customize the modeling, solving and cutting process.

## Key Concepts

`AbstractData`: The data structure that stores the problem data. For large-scale problems, it is recommended to use the `AbstractData` interface to store the problem data. Then the user can use the data to build the model.

`SolutionProcedure`: Itertively or callback function to process the solution.

`CutStrategy`: The strategy to generate the cuts, such as classical Benders cut, unified Benders cut, etc.

`AbstractMasterProblem`: The abstract type for the master problem.

`AbstractSubProblem`: The abstract type for the subproblem.

## Modeling

Unlike algebraic algorithms, Benders decomposition is typically used for large-scale problems. In this approach, researchers do not rely on abstract matrix representations to compute key elements in the solving process.
Instead, researchers usually explicitly formulate both the master problem and the subproblem. The partitioning of variables and constraints can vary significantly depending on the application.

Then the next question is how to build the model. Which information should we use to build the model? Our answer is `AbstractData` and `CutStrategy`. For different applications, the `AbstractData` is different. For the same application, the `CutStrategy` can be different.
We only need to provide the `AbstractData` and `CutStrategy` to build the model.
So we have two functions to build the model:
- Use `create_master_problem(data::AbstractData, cut_strategy::CutStrategy)` to build `master_problem <: AbstractMasterProblem`.
- Use `create_sub_problem(data::AbstractData, cut_strategy::CutStrategy)` to build `sub_problem <: AbstractSubProblem`.

We provide a general standard structure for master problem and subproblem for different applications. General structures:
```julia
struct GeneralMasterProblem <: AbstractMasterProblem
    model::Model # JuMP model
    var::Dict # variable dictionary
    obj_value::Float64 # objective value
    x_value::Vector{Float64} # integer variables
    t_value::Union{Float64,Vector{Float64}} # continuous variables
end
```
and
```julia
struct GeneralSubProblem <: AbstractSubProblem
    model::Model # JuMP model
    fixed_x_constraints::Vector{ConstraintRef} # fixed integer variables constraints
    other_constraints::Vector{ConstraintRef} # other constraints
end
```
But if user want to use their own unique model, they can define their own structures and multiple dispatch the `create_master_problem` and `create_sub_problem` function.

For example, the fat knapsack cut method proposed by Fischetti in uncapacitated facility location problem, we don't need to build a traditional subproblem. We define a new structure storing necessary information to represent the subproblem:
```julia
mutable struct KnapsackUFLPSubProblem <: AbstractUFLPSubProblem
    sorted_cost_demands::Vector{Vector{Float64}}
    sorted_indices::Vector{Vector{Int}}
    selected_k::Dict
end
```

If you want to use the Disjunctive Benders Cuts we proposed in our paper, you can use `create_dcglp(data::AbstractData, cut_strategy::CutStrategy)` to build `dcglp <: AbstractDCGLP`.

## Environment

Here we introduce a new concept: `Environment`, borrowed from Reinforcement Learning.
The `BendersEnv` is a structure that stores the modeling information:
```julia
mutable struct BendersEnv
    master::AbstractMasterProblem
    sub::Union{AbstractSubProblem, Vector{AbstractSubProblem}} # for deterministic problems, sub is a single subproblem; for stochastic problems, sub is a vector of subproblems
    dcglp::Union{Nothing, DCGLP}  # Optional component
end
```

## Solution Procedure
In textbooks, Benders decomposition is typically presented as an iterative algorithm.
However, in research applications, the Branch-and-Cut method is often used to obtain the final results. In this approach, researchers leverage callback methods to incorporate lazy constraints, user cuts and heuristic methods, providing flexibility in implementation. This allows for various strategies to enhance the efficiency of the decomposition process.

In this package, we provide `Sequential()`, `Callback()`, `StochasticSequential()`, `StochasticCallback()` to solve the problem.

The corresponding function is `solve!(env::BendersEnv, solution_procedure::SolutionProcedure, params::Params)`.

For example, the iterative solution procedure is:
```julia
function solve!(env::BendersEnv, ::Sequential, cut_strategy::CutStrategy, params::BendersParams)
    log = BendersIterationLog()
    state = BendersState()

    while true
        state.iteration += 1

        # Solve master problem
        master_time = @elapsed begin
            solve_master!(env.master)
            state.LB = env.master.obj_value
        end
        log.master_time += master_time

        # Solve sub problem
        sub_time = @elapsed begin
            solve_sub!(env.sub, env.master.x_value)
            cuts, sub_obj_val = generate_cuts(env, cut_strategy)
            update_upper_bound_and_gap!(state, env, sub_obj_val)
        end
        log.sub_time += sub_time

        # Update state and record information
        record_iteration!(log, state)

        params.verbose && print_iteration_info(state, log)

        # Check termination criteria
        is_terminated(state, params, log) && break

        # Generate and add cuts
        for cut in cuts
            @constraint(env.master.model, 0 >= cut)
        end
    end
end
```

The callback method:
```julia
function solve!(env::BendersEnv, ::Callback, cut_strategy::CutStrategy, params::BendersParams)

    start_time = time()
    function lazy_callback(cb_data)
        status = JuMP.callback_node_status(cb_data, env.master.model)
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER
            env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
            env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])

            solve_sub!(env.sub, env.master.x_value)
            cuts, sub_obj_value = generate_cuts(env, cut_strategy)
            add_cuts!(env, cuts, sub_obj_value, cb_data)
        end
    end

    function user_cut_callback(cb_data)
        ... # user cut callback
    end

    # Use the closure callbacks
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback)
    set_attribute(env.master.model, MOI.UserCutCallback(), user_cut_callback)
    MOI.set(env.master.model, MOI.RelativeGapTolerance(), params.gap_tolerance)
    set_time_limit_sec(env.master.model, params.time_limit)
    set_optimizer_attribute(env.master.model, MOI.Silent(), false)
    JuMP.optimize!(env.master.model)
end
```

The user can customize the solution procedure by defining a new `SolutionProcedure` and multiple dispatch the `solve!` function.

## Cut Strategy

The cut strategy is the key to the Benders decomposition algorithm. The cut strategy determines how to generate the cuts.

We provide some general cut strategies for different applications, such as classical cut, unified cut, disjunctive cut, etc.

The key function is `generate_cut(env::BendersEnv, cut_strategy::CutStrategy)`, which generates the cuts and add them to the model.

For example, the classical cut is:
```julia
function generate_cuts(env::BendersEnv, ::ClassicalCut)
    (coefficients_t, coefficients_x, constant_term), sub_obj_val = generate_cut_coefficients(env.sub, env.master.x_value, ClassicalCut())

    cut = @expression(env.master.model,
        constant_term + dot(coefficients_x, env.master.var[:x]) + coefficients_t * env.master.var[:t])

    return cut, sub_obj_val
end
```

If you want to use the fat knapsack cut:
```julia
function generate_cuts(env::BendersEnv, cut_strategy::Union{FatKnapsackCut, SlimKnapsackCut})

    # return vector
    critical_pairs, obj_values = generate_cut_coefficients(env.sub, env.master.x_value, cut_strategy)

    cuts = Vector{Any}(undef, length(critical_pairs))
    for (index, critical_item) in critical_pairs
        cuts[index] = build_cut(env.master, env.sub, (index, critical_item), cut_strategy)
    end

    return cuts, obj_values
end
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

