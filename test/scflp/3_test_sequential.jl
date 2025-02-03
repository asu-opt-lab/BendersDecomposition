using Test
using JuMP
using CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition

# global_logger(ConsoleLogger(stderr, Logging.Warn))

@testset "SCFLP Sequential Benders Tests" begin
    solver = "CPLEX"
    # solver = :Gurobi
    data = read_stochastic_capacited_facility_location_problem("f25-c50-s64-r10-1")

    # Create and solve MIP reference model
    milp = create_milp(data)
    set_optimizer(milp.model, CPLEX.Optimizer)
    optimize!(milp.model)
    mip_obj = objective_value(milp.model)

    # Test different cut strategies
    loop_strategy = StochasticSequential()
    cut_strategies = Dict(
        "Standard" => ClassicalCut(),
        "Knapsack" => KnapsackCut()
    )
    params = BendersParams(
        600.0,
        0.00001,
        solver,
        Dict("solver" => solver),
        Dict("solver" => solver),
        Dict(),
        true
    )

    benders_UB = Dict()
    benders_LB = Dict()
    
    for (name, strategy) in cut_strategies
        result = run_Benders(data, loop_strategy, strategy, params)
        benders_LB[name] = result[end, :LB]
        benders_UB[name] = result[end, :UB]
        
        @test isapprox(mip_obj, benders_LB[name], atol=0.1)
        @test isapprox(mip_obj, benders_UB[name], atol=0.1)
    end

    # Print results
    @printf("MIP: %.4f | ", mip_obj)
    for (name, lb) in benders_LB
        @printf("%s_LB: %.4f | %s_UB: %.4f | ", name, lb, name, benders_UB[name])
    end
    println()
end
