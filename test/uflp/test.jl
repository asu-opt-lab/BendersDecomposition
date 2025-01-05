using Test
using JuMP
using Gurobi, CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition


solver = "CPLEX"
    
# instances=[
#     "ga500a-1" "ga500a-2" "ga500a-3" "ga500a-4" "ga500a-5"
#     "ga500b-1" "ga500b-2" "ga500b-3" "ga500b-4" "ga500b-5"
#     "ga750a-1" "ga750a-2" "ga750a-3" "ga750a-4" "ga750a-5"
#     "ga750b-1" "ga750b-2" "ga750b-3" "ga750b-4" "ga750b-5"
#     "ga750c-1" "ga750c-2" "ga750c-3" "ga750c-4" "ga750c-5"
# ]

instances=[
    "ga500a-1"
]


for i in instances
    # data = read_Simple_data("ga250a-3")
    data = read_Simple_data("$i")
    @info i

    disjunctive_system = DisjunctiveCut(FatKnapsackCut(), L1Norm(), PureDisjunctiveCut(), true, true, false,false)
    @info disjunctive_system
    params = BendersParams(
        7200.0,
        1e-5, # *100 already
        solver,
        Dict("solver" => solver),
        Dict("solver" => solver),
        Dict("solver" => solver),
        # Dict(:solver => :Gurobi),
        true
        # false
    )

    # result = run_Benders(data, Sequential(), ClassicalCut(), params)
    result = run_Benders(data, Sequential(), FatKnapsackCut(), params)
    # result = run_Benders(data, Callback(), FatKnapsackCut(), params)
    # result = run_Benders(data, Sequential(), disjunctive_system, params)
    # result = run_Benders(data, Callback(), disjunctive_system, params)
end