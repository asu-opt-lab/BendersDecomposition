# Main Function

This package implements Benders decomposition.

## Core Function
`run_Benders(data::AbstractData, loop_strategy::SolutionProcedure, cut_strategy::CutStrategy, params::BendersParams)` 
is the main entry point for running the Benders decomposition algorithm. It requires four arguments:

```julia
function run_Benders(data::AbstractData, loop_strategy::SolutionProcedure, cut_strategy::CutStrategy, params::BendersParams)
    env = BendersEnv(data, cut_strategy, params)
    solve!(env, loop_strategy, cut_strategy, params)
end
```

## Environment Structure

The `BendersEnv` type contains all components needed for Benders decomposition:

```julia
mutable struct BendersEnv
    data::AbstractData                                        # Problem instance data
    master::AbstractMasterProblem                            # Master problem
    sub::Union{AbstractSubProblem, Vector{AbstractSubProblem}} # Subproblem(s)
    dcglp::Union{Nothing, DCGLP}                            # Optional DCGLP component
end
```

### Environment Creation

The environment is initialized using:

```julia
function BendersEnv(data::AbstractData, cut_strategy::CutStrategy, params::BendersParams)
    # Create and configure master problem
    master = create_master_problem(data, cut_strategy)
    assign_attributes!(master.model, params.master_attributes)
    
    # Create and configure subproblem(s)
    sub = create_sub_problem(data, cut_strategy)
    if sub isa AbstractSCFLPSubProblem
        for scenario_sub in sub.sub_problems
            assign_attributes!(scenario_sub.model, params.sub_attributes)
        end
    else
        assign_attributes!(sub.model, params.sub_attributes)
    end
    
    # Create DCGLP if using disjunctive cuts
    dcglp = if cut_strategy isa DisjunctiveCut
        dcglp = create_dcglp(data, cut_strategy)
        assign_attributes!(dcglp.model, params.dcglp_attributes)
        dcglp
    else
        nothing
    end
    
    return BendersEnv(data, master, sub, dcglp)
end
```

### Components

1. **Master Problem**: Created via `create_master_problem(data, cut_strategy)`. The structure depends on the chosen cut strategy.

2. **Subproblem**: Created via `create_sub_problem(data, cut_strategy)`. Can be either a single problem or a vector of problems for scenario-based formulations.

3. **DCGLP** (Optional): Created via `create_dcglp(data, cut_strategy)` when using disjunctive cuts. Structure:

```julia
mutable struct DCGLP <: AbstractDCGLP
    model::Model
    γ_constraints::Dict{Symbol,Any}
    γ_values::Vector{Tuple{Float64, Vector{Float64}, Union{Float64, Vector{Float64}}}}
    disjunctive_inequalities_constraints::Vector{ConstraintRef}
    dcglp_constraints::Any
    master_cuts::Any
end
```

The DCGLP is created using:

```julia
function create_dcglp(data::AbstractData, disjunction_system::DisjunctiveCut)
    model = create_base_model(data, disjunction_system.norm_type)
    add_problem_specific_constraints!(model, data, disjunction_system.norm_type)
    add_t_constraints!(model, data, disjunction_system.base_cut_strategy, disjunction_system.norm_type)
    add_norm_specific_components!(model, data, disjunction_system.base_cut_strategy, disjunction_system.norm_type)

    return DCGLP(model, collect_constraints(model), [], [], [], [])
end
```

# solve!

## Loop Strategy

```julia
function solve!(env::BendersEnv, loop::SolutionProcedure, cut_strategy::CutStrategy, params::BendersParams)
    error("solve! not implemented for $(typeof(loop)) with $(typeof(cut_strategy))")
end
```

For different loop strategy, the `solve!` function is different. We have three loop strategies:

1. Sequential
2. Callback
3. StochasticSequential

The `solve!` function for each loop strategy is implemented in the corresponding file and only related to the loop strategy. 

## Cut Strategy

```julia
function generate_cuts(env::BendersEnv, cut_strategy::CutStrategy)
    error("generate_cuts not implemented for strategy type $(typeof(cut_strategy))")
end
```

Because in dcglp, we need some information on dual variables, we have the following function:

```julia
function generate_cut_coefficients(sub::AbstractSubProblem, x_value::Vector{Float64}, cut_strategy::CutStrategy)
    error("generate_cut_coefficients not implemented for subproblem type $(typeof(sub)) and strategy $(typeof(cut_strategy))")
end
```
For disjunctive cut, it's little different. 

```julia
function generate_cuts(env::BendersEnv, cut_strategy::DisjunctiveCut)

    sub_obj_val = get_subproblem_value(env) 

    disjunctive_inequality = select_disjunctive_inequality(env.master.x_value)

    update_dcglp!(env.dcglp, disjunctive_inequality, cut_strategy)
    
    solve_dcglp!(env, cut_strategy)
    
    cuts = merge_cuts(env, cut_strategy)

    return cuts, sub_obj_val
end
```

and

```julia
function update_dcglp!(dcglp::DCGLP, disjunctive_inequality::Tuple{Vector{Int}, Int}, disjunction_system::DisjunctiveCut)
    replace_disjunctive_inequality!(dcglp, disjunctive_inequality, disjunction_system.norm_type)
    update_added_benders_constraints!(dcglp, disjunction_system)
    add_disjunctive_cut!(dcglp)
end 
```

## Solution Process

The solution process consists of two main components: the loop strategy and the cut generation strategy.

### Loop Strategies

The `solve!` function implements different solution procedures based on the chosen loop strategy:

```julia
function solve!(env::BendersEnv, loop::SolutionProcedure, cut_strategy::CutStrategy, params::BendersParams)
    error("solve! not implemented for $(typeof(loop)) with $(typeof(cut_strategy))")
end
```

Three main loop strategies are available:
1. **Sequential**: Traditional iterative approach
2. **Callback**: Solver-based lazy constraint generation
3. **StochasticSequential**: Specialized for stochastic problems

Each strategy is implemented in its corresponding file with strategy-specific logic.

### Cut Generation

Cut generation is handled through two main interfaces:

1. **Cut Generation**:
```julia
function generate_cuts(env::BendersEnv, cut_strategy::CutStrategy)
    error("generate_cuts not implemented for strategy type $(typeof(cut_strategy))")
end
```

2. **Cut Coefficients** (for dual information):
```julia
function generate_cut_coefficients(sub::AbstractSubProblem, x_value::Vector{Float64}, cut_strategy::CutStrategy)
    error("generate_cut_coefficients not implemented for subproblem type $(typeof(sub)) and strategy $(typeof(cut_strategy))")
end
```

#### Disjunctive Cuts

For disjunctive cuts, the generation process involves the DCGLP:

```julia
function generate_cuts(env::BendersEnv, cut_strategy::DisjunctiveCut)
    # Get subproblem objective value
    sub_obj_val = get_subproblem_value(env) 

    # Select appropriate disjunctive inequality
    disjunctive_inequality = select_disjunctive_inequality(env.master.x_value)

    # Update and solve DCGLP
    update_dcglp!(env.dcglp, disjunctive_inequality, cut_strategy)
    solve_dcglp!(env, cut_strategy)
    
    # Combine results into final cuts
    cuts = merge_cuts(env, cut_strategy)

    return cuts, sub_obj_val
end
```

The DCGLP update process:
```julia
function update_dcglp!(dcglp::DCGLP, disjunctive_inequality::Tuple{Vector{Int}, Int}, disjunction_system::DisjunctiveCut)
    replace_disjunctive_inequality!(dcglp, disjunctive_inequality, disjunction_system.norm_type)
    update_added_benders_constraints!(dcglp, disjunction_system)
    add_disjunctive_cut!(dcglp)
end 
```

