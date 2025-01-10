export StandardSCFLPSubProblem

abstract type AbstractSCFLPSubProblem <: AbstractSubProblem end

mutable struct StandardSCFLPSubProblem <: AbstractSCFLPSubProblem
    sub_problems::Vector{StandardCFLPSubProblem}
end

mutable struct KnapsackSCFLPSubProblem <: AbstractSCFLPSubProblem
    sub_problems::Vector{KnapsackCFLPSubProblem}
end

# Specialized create_sub_problem functions
function create_sub_problem(data::SCFLPData, ::ClassicalCut)
    @debug "Building Subproblem for SCFLP (Standard)"
    sub_problems = Vector{StandardCFLPSubProblem}()
    for scenario in 1:data.n_scenarios
        data_scenario = CFLPData(data.n_facilities, data.n_customers, data.capacities, data.demands[scenario], data.fixed_costs, data.costs)
        model, fixed_x, other = _create_base_sub_problem(data_scenario)
        push!(sub_problems, StandardCFLPSubProblem(model, fixed_x, other))
    end
    return StandardSCFLPSubProblem(sub_problems)
end

function create_sub_problem(data::SCFLPData, ::KnapsackCut)
    @debug "Building Subproblem for SCFLP (Knapsack)"
    sub_problems = Vector{KnapsackCFLPSubProblem}()
    for scenario in 1:data.n_scenarios
        data_scenario = CFLPData(data.n_facilities, data.n_customers, data.capacities, data.demands[scenario], data.fixed_costs, data.costs)
        model, fixed_x, other, demand = _create_base_sub_problem(data_scenario)
        facility_knapsack_info = FacilityKnapsackInfo(data_scenario.costs, data_scenario.demands, data_scenario.capacities)
        push!(sub_problems, KnapsackCFLPSubProblem(model, fixed_x, other, demand, facility_knapsack_info))
    end
    return KnapsackSCFLPSubProblem(sub_problems)
end

