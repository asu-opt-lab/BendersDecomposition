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
    # solver = :Gurobi
    # instances = [1:66; 68:71]
    instances = [20]
    for i in instances
        @testset "Instance: p$i" begin
            # Load data
            # data = read_cflp_benchmark_data("p$i")
            data = read_GK_data("f100-c100-r5-1")
            # Create and solve MIP reference model
            milp = create_milp(data)
            set_optimizer(milp.model, CPLEX.Optimizer)
            optimize!(milp.model)
            mip_obj = objective_value(milp.model)
            
            # Test different cut strategies
            loop_strategy = Callback()
            cut_strategies = Dict(
                # "Classical" => ClassicalCut(),
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
            # benders_UB = Dict()
            # benders_LB = Dict()
            for (name, strategy) in cut_strategies
                obj_value, _ = run_Benders(data, loop_strategy, strategy, params)
                # @test isapprox(mip_obj, obj_value, atol=0.1)
            end
            # @printf("Instance: p%d | %-5s: %.2f | ", i, "MIP", mip_obj)
            # for (cut_type, lb) in benders_LB
            #     @printf("%-5s: %.2f | ", string(cut_type), lb)
            # end
            # println()
        end
    end
end
