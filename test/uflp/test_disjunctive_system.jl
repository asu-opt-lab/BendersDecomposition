using Test
using JuMP
using Gurobi, CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition

global_logger(ConsoleLogger(stderr, Logging.Warn))

@testset "UFLP Disjunctive System Tests" begin
    solver = "CPLEX"
    
    # Test on a few representative instances
    for i in [1:66;68:71]
        @testset "Instance: p$(i)" begin
            # Load UFLP data
            data = read_uflp_benchmark_data("p$(i)")
            
            # Solve using standard MIP model for comparison
            milp = create_milp(data)
            set_optimizer(milp.model, CPLEX.Optimizer)
            optimize!(milp.model)
            mip_objective = objective_value(milp.model)
            
            # disjunctive_system = DisjunctiveCut(FatKnapsackCut(), L1Norm(), PureDisjunctiveCut(), true, true, true,true)
            # disjunctive_system = DisjunctiveCut(SlimKnapsackCut(), L1Norm(), PureDisjunctiveCut(), true, true, true,true)
            disjunctive_system = DisjunctiveCut(ClassicalCut(), L1Norm(), PureDisjunctiveCut(), true, true, true,true)
            loop_strategy = Sequential()
            
            params = BendersParams(7200.0, 1e-9, solver, Dict("solver" => solver), Dict("solver" => solver), Dict("solver" => solver), true) 
            
            result = run_Benders(data, loop_strategy, disjunctive_system, params)
            disjunctive_objective = result[end, :UB]
            
            # Compare results
            @test isapprox(mip_objective, disjunctive_objective, rtol=1e-4)
            
            # Print results
            @printf("Instance: p%d | MIP: %.4f | Disjunctive: %.4f | Iterations: %d | Time: %.2f\n",
                    i, mip_objective, disjunctive_objective, result[end, :iter], result[end, :total_time])
        end
    end
end
