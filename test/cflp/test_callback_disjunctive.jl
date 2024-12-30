using Test
using JuMP
using CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition

# global_logger(ConsoleLogger(stderr, Logging.Warn))

@testset "CFLP Callback Benders Tests" begin
    solver = "Gurobi"
    # solver = :CPLEX
    instances = [1:66; 68:71]
    # instances = [2]
    for i in instances
        @testset "Instance: p$i" begin
            # Load data
            data = read_cflp_benchmark_data("p$i")
            
            # Create and solve MIP reference model
            milp = create_milp(data)
            set_optimizer(milp.model, CPLEX.Optimizer)
            optimize!(milp.model)
            mip_obj = objective_value(milp.model)
            
            # Test different cut strategies
            loop_strategy = Callback()
            disjunctive_system = DisjunctiveCut(KnapsackCut(), LInfNorm(), PureDisjunctiveCut(), true, false,false,false)
            # disjunctive_system = DisjunctiveCut(KnapsackCut(), LInfNorm(), PureDisjunctionCut(), true, false,false,false)
            
            params = BendersParams(
                600.0,
                0.000001,
                solver,
                Dict("solver" => solver),
                Dict("solver" => solver),
                Dict("solver" => solver),
                # Dict(:solver => :Gurobi),
                # true
                false
            )
            
            
            
            result1, result2 = run_Benders(data, loop_strategy, disjunctive_system, params)
            disjunctive_LB = result2[end, :objective_bound] # Use LB instead of UB
            disjunctive_UB = result2[end, :objective_value]
            # Compare results
            @test isapprox(mip_obj, disjunctive_LB, atol=1e-1)
            # Test if LB is approximately equal to UB
            @test isapprox(disjunctive_LB, disjunctive_UB, atol=1e-1)
            # Print results
            @printf("Instance: p%d | MIP: %.4f | Disjunctive_LB: %.4f | Disjunctive_UB: %.4f |Iterations: %d | Time: %.2f\n",
                    i, mip_obj, disjunctive_LB, disjunctive_UB, result2[end, :node_count], result2[end, :elapsed_time])
        end
    end
end
