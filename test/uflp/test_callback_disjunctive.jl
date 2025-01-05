using Test
using JuMP
using Gurobi, CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition


# solver = "CPLEX"
    
# # data = read_uflp_benchmark_data("p1")
# data = read_Simple_data("ga250a-3")

# disjunctive_system = DisjunctiveCut(FatKnapsackCut(), LInfNorm(), PureDisjunctiveCut(), true, true, true,false)

# params = BendersParams(
#     7200.0,
#     1e-5, # *100 already
#     solver,
#     Dict("solver" => solver),
#     Dict("solver" => solver),
#     Dict("solver" => solver),
#     # Dict(:solver => :Gurobi),
#     true
#     # false
# )

# # result = run_Benders(data, Sequential(), FatKnapsackCut(), params)
# # result = run_Benders(data, Callback(), FatKnapsackCut(), params)
# result = run_Benders(data, Sequential(), disjunctive_system, params)
# # result = run_Benders(data, Callback(), disjunctive_system, params)

"number of cuts added to master problem
use_two_sided_cuts / include_master_cuts / reuse_dcglp
1) true/true/true: 50, 7 iterations
2) true/false/false: at least add j (# of subproblems) cuts, many iterations
3) true/true/false: 77, 7 iterations
"


global_logger(ConsoleLogger(stderr, Logging.Warn))

@testset "UFLP Disjunctive System Tests" begin
    solver = "CPLEX"

    instances = [1:66; 68:71]

    # Test on a few representative instances
    for i in instances
        @testset "Instance: p$(i)" begin
            # Load UFLP data
            data = read_uflp_benchmark_data("p$(i)")

            # Solve using standard MIP model for comparison
            milp = create_milp(data)
            set_optimizer(milp.model, CPLEX.Optimizer)
            optimize!(milp.model)
            mip_objective = objective_value(milp.model)
            
            disjunctive_system = DisjunctiveCut(FatKnapsackCut(), LInfNorm(), PureDisjunctiveCut(), true, true, true,false)
            loop_strategy = Callback()
            
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
            result1, result2 = run_Benders(data, loop_strategy, disjunctive_system, params)
            disjunctive_LB = result2[end, :objective_bound] # Use LB instead of UB
            disjunctive_UB = result2[end, :objective_value]
            # Compare results
            @test isapprox(mip_objective, disjunctive_LB, atol=1e-1)
            # Test if LB is approximately equal to UB
            @test isapprox(disjunctive_LB, disjunctive_UB, atol=1e-1)
            # Print results
            @printf("Instance: p%d | MIP: %.4f | Disjunctive_LB: %.4f | Disjunctive_UB: %.4f |Iterations: %d | Time: %.2f\n",
                    i, mip_objective, disjunctive_LB, disjunctive_UB, result2[end, :node_count], result2[end, :elapsed_time])
        end
    end
end