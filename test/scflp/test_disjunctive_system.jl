using Test
using JuMP
using CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition

# global_logger(ConsoleLogger(stderr, Logging.Warn))

@testset "SCFLP Sequential Benders Tests" begin
    # solver = "CPLEX"
    solver = "Gurobi"
    data = read_stochastic_capacited_facility_location_problem("f50-c100-s256-r10-1")

    # Create and solve MIP reference model
    # milp = create_milp(data)
    # set_optimizer(milp.model, CPLEX.Optimizer)
    # optimize!(milp.model)
    # mip_obj = objective_value(milp.model)

    # Test different cut strategies
    loop_strategy = StochasticSequential()
    disjunctive_system = DisjunctiveCut(KnapsackCut(), LInfNorm(), PureDisjunctiveCut(), true, true, false,true)
    # disjunctive_system = DisjunctiveCut(ClassicalCut(), L2Norm(), PureDisjunctiveCut(), true, true, false,true)
    params = BendersParams(
        600.0,
        0.00001,
        solver,
        Dict("solver" => solver),
        Dict("solver" => solver),
        Dict("solver" => solver),
        true
    )

    
    result = run_Benders(data, loop_strategy, disjunctive_system, params)
    disjunctive_LB = result[end, :LB] 
    disjunctive_UB = result[end, :UB]
    
    @test isapprox(mip_obj, disjunctive_LB, atol=0.1)
    @test isapprox(mip_obj, disjunctive_UB, atol=0.1)

    # Print results
    @printf("MIP: %.4f | Disjunctive_LB: %.4f | Disjunctive_UB: %.4f |Iterations: %d | Time: %.2f\n",
    mip_obj, disjunctive_LB, disjunctive_UB, result[end, :iter], result[end, :total_time])
    println()
end
