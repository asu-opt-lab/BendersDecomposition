using Test
using JuMP
using CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition


# global_logger(ConsoleLogger(stderr, Logging.Warn))



@testset "CFLP Sequential Benders Tests" begin
    solver = "Gurobi"
    # solver = :Gurobi
    K = "04"
    instances = ["r$K.$j.dow" for j in 1:9]
    for i in instances
        @testset "Instance: $i" begin
            # Load data
            data = read_mcndp_instance(i)
            
            # Create and solve MIP reference model
            milp = create_milp(data)
            set_optimizer(milp.model, CPLEX.Optimizer)
            optimize!(milp.model)
            mip_obj = objective_value(milp.model)
            @info mip_obj
            # Test different cut strategies
            loop_strategy = Callback()
            cut_strategies = Dict(
                "Standard" => ClassicalCut(),
                "Knapsack" => KnapsackCut()
            )
            params = BendersParams(
                600.0,
                0.00001,
                solver,
                Dict("solver" => solver),
                Dict("solver" => solver, "InfUnbdInfo" => 1),
                Dict("solver" => solver),
                true
            )
            benders_UB = Dict()
            benders_LB = Dict()
            for (name, strategy) in cut_strategies
                obj_value, _ = run_Benders(data, loop_strategy, strategy, params)
                @test isapprox(mip_obj, obj_value, atol=0.1)
            end
            @printf("Instance: %s | MIP: %.2f | ", i, mip_obj)
            for (cut_type, lb) in benders_LB
                @printf("%s: %.2f | ", string(cut_type), lb)
            end
            println()
        end
    end
end
