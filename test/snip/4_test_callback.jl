using Test
using JuMP
using CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition


# global_logger(ConsoleLogger(stderr, Logging.Warn))



@testset "SNIP Sequential Benders Tests" begin
    solver = "Gurobi"
    # solver = :Gurobi

    instances = [0, 1, 2, 3, 4]
    # instances = [0]
    for i in instances
        @testset "Instance: $i" begin
            # Load data
            data = read_snip_data(i, 1, 30.0)
            
            # Create and solve MIP reference model
            milp = create_milp(data)
            set_optimizer(milp.model, CPLEX.Optimizer)
            optimize!(milp.model)
            mip_obj = objective_value(milp.model)
            @info mip_obj
            
            # Test different cut strategies
            # loop_strategy = StochasticSequential()
            loop_strategy = StochasticCallback()
            cut_strategies = Dict(
                "Classical" => ClassicalCut()
            )
            params = BendersParams(
                1000.0,
                0.001,
                solver,
                Dict("solver" => solver),
                Dict("solver" => solver),
                # Dict("solver" => solver, "InfUnbdInfo" => 1),
                Dict("solver" => solver),
                true
            )
            benders_UB = Dict()
            benders_LB = Dict()
            for (name, strategy) in cut_strategies
                result = run_Benders(data, loop_strategy, strategy, params)
                # @info mip_obj
                benders_LB[name] = result[end, :LB]
                benders_UB[name] = result[end, :UB]
                @test isapprox(mip_obj, result[end, :LB], atol=0.1)
                @test isapprox(mip_obj, result[end, :UB], atol=0.1)
            end
            @printf("Instance: %s | MIP: %.2f | ", i, mip_obj)
            for (cut_type, lb) in benders_LB
                @printf("%s: %.2f | ", string(cut_type), lb)
            end
            println()
        end
    end
end
