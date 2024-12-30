using Test
using JuMP
using Gurobi, CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition

# global_logger(ConsoleLogger(stderr, Logging.Warn))

@testset "CFLP Disjunctive System Tests" begin
    # solver = :Gurobi
    solver = "CPLEX"
    
    # Test on a few representative instances
    # for i in [1:66;68:71]
    # for i in 29:66
    for i in [14,30,60,71]
        @testset "Instance: p$(i)" begin
            # Load CFLP data
            data = read_cflp_benchmark_data("p$(i)")
            
            # Solve using standard MIP model for comparison
            milp = create_milp(data)
            set_optimizer(milp.model, CPLEX.Optimizer)
            optimize!(milp.model)
            mip_objective = objective_value(milp.model)
            

            # disjunctive_system = DisjunctiveCut(ClassicalCut(), L1Norm(), PureDisjunctiveCut(), true, true, false,true)
            disjunctive_system = DisjunctiveCut(KnapsackCut(), LInfNorm(), PureDisjunctiveCut(), true, false,false,false)
            
            params = BendersParams(
                60.0,
                1e-5, # *100 already
                solver,
                Dict("solver" => solver),
                Dict("solver" => solver),
                Dict("solver" => solver),
                # Dict(:solver => :Gurobi),
                true
                # false
            )
            result = run_Benders(data, Sequential(), disjunctive_system, params)
            disjunctive_LB = result[end, :LB] # Use LB instead of UB
            disjunctive_UB = result[end, :UB]
            # Compare results
            @test isapprox(mip_objective, disjunctive_LB, atol=1e-1)
            # Test if LB is approximately equal to UB
            @test isapprox(disjunctive_LB, disjunctive_UB, atol=1e-1)
            # Print results
            @printf("Instance: p%d | MIP: %.4f | Disjunctive_LB: %.4f | Disjunctive_UB: %.4f |Iterations: %d | Time: %.2f\n",
                    i, mip_objective, disjunctive_LB, disjunctive_UB, result[end, :iter], result[end, :total_time])
        end
    end
end
