using Test
using JuMP
using Gurobi, CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition

global_logger(ConsoleLogger(stderr, Logging.Warn))

@testset "UFLP Disjunctive System Tests" begin
    solver = :CPLEX
    
    # Test on a few representative instances
    for i in [1:66;68:71]
        @testset "Instance: p$(i)" begin
            # Load UFLP data
            data = read_uflp_benchmark_data("p$(i)")
            
            # Solve using standard MIP model for comparison
            milp = create_milp(data)
            assign_solver!(milp.model, solver)
            optimize!(milp.model)
            mip_objective = objective_value(milp.model)
            
            disjunctive_system = DisjunctionSystem(L1Norm(), FatKnapsackCut())
            # Solve using disjunctive Benders
            master = create_master_problem(data, disjunctive_system._cut_strategy)
            relax_integrality(master.model)
            assign_solver!(master.model, solver)
            sub = create_sub_problem(data, disjunctive_system._cut_strategy)
            # assign_solver!(sub.model, solver)
            
            params = BendersParams(600.0, 0.001)  # 10 minutes time limit, 0.1% gap tolerance
            algo = SequentialBenders(data, master, sub, disjunctive_system, params)
            
            result = solve!(algo)
            disjunctive_objective = result[end, :UB]
            
            # Compare results
            @test isapprox(mip_objective, disjunctive_objective, rtol=1e-4)
            
            # Print results
            @printf("Instance: p%d | MIP: %.4f | Disjunctive: %.4f | Iterations: %d | Time: %.2f\n",
                    i, mip_objective, disjunctive_objective, result[end, :iter], result[end, :elapsed_time])
        end
    end
end
