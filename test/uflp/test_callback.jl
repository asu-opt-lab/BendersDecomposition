using Test
using JuMP
using CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition


# global_logger(ConsoleLogger(stderr, Logging.Warn))



@testset "CFLP Sequential Benders Tests" begin
    solver = "CPLEX"
    # instances = [1:66; 68:71]
    instances = [2]
    for i in instances
        @testset "Instance: p$i" begin
            # Load data
            data = read_uflp_benchmark_data("p$i")
            
            # Create and solve MIP reference model
            milp = create_milp(data)
            set_optimizer(milp.model, CPLEX.Optimizer)
            optimize!(milp.model)
            mip_obj = objective_value(milp.model)
            
            # Test different cut strategies
            loop_strategy = Callback()
            cut_strategies = Dict(
                "Standard" => ClassicalCut(),
                # "Knapsack" => FatKnapsackCut()
            )
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
            benders_UB = Dict()
            benders_LB = Dict()
            for (name, strategy) in cut_strategies
                result1, result2 = run_Benders(data, loop_strategy, strategy, params)

                benders_LB[name] = result2[end, :objective_bound]
                benders_UB[name] = result2[end, :objective_value]
                @test isapprox(mip_obj, result2[end, :objective_bound], atol=0.1)
                @test isapprox(mip_obj, result2[end, :objective_value], atol=0.1)
            end
            @printf("Instance: p%d | %-5s: %.2f | ", i, "MIP", mip_obj)
            for (cut_type, lb) in benders_LB
                @printf("%-5s: %.2f | ", string(cut_type), lb)
            end
            println()
        end
    end
end
