export StandardSCFLPSubProblem, KnapsackSCFLPSubProblem


abstract type AbstractSCFLPSubProblem <: AbstractSubProblem end



"""
    StandardSCFLPSubProblem <: AbstractSCFLPSubProblem

A mutable struct representing the standard subproblem formulation for the Stochastic Capacitated Facility Location Problem (SCFLP).

# Fields
- `sub_problems::Vector{StandardCFLPSubProblem}`: Vector of scenario-specific CFLP subproblems

# Related Functions
    create_sub_problem(data::SCFLPData, ::ClassicalCut)
"""
mutable struct StandardSCFLPSubProblem <: AbstractSCFLPSubProblem
    sub_problems::Vector{StandardCFLPSubProblem}
end

"""
    KnapsackSCFLPSubProblem <: AbstractSCFLPSubProblem

A mutable struct representing the knapsack-based subproblem formulation for the SCFLP.

# Fields
- `sub_problems::Vector{KnapsackCFLPSubProblem}`: Vector of scenario-specific knapsack CFLP subproblems

# Related Functions
    create_sub_problem(data::SCFLPData, ::KnapsackCut)
"""
mutable struct KnapsackSCFLPSubProblem <: AbstractSCFLPSubProblem
    sub_problems::Vector{KnapsackCFLPSubProblem}
end

function create_scenario_data(data::SCFLPData, scenario::Int)
    return CFLPData(
        data.n_facilities,
        data.n_customers,
        data.capacities,
        data.demands[scenario],
        data.fixed_costs,
        data.costs
    )
end

function create_sub_problem(data::SCFLPData, cut_strategy::Union{ClassicalCut, KnapsackCut})
    @debug "Building Subproblem for SCFLP (Standard)"
    sub_problems = map(1:data.n_scenarios) do scenario
        data_scenario = create_scenario_data(data, scenario)
        create_sub_problem(data_scenario, cut_strategy)
    end
    return _wrap_sub_problems(sub_problems)
end

_wrap_sub_problems(sub_problems::Vector{StandardCFLPSubProblem}) = StandardSCFLPSubProblem(sub_problems)
_wrap_sub_problems(sub_problems::Vector{KnapsackCFLPSubProblem}) = KnapsackSCFLPSubProblem(sub_problems)

