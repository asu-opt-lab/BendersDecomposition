using Test
using JuMP
using Gurobi, CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition


global_logger(ConsoleLogger(stderr, Logging.Warn))



@testset "UFLP Instances Comparison" begin
    # Get all UFLP instance files
    solver = "CPLEX"
    instances = [1:66; 68:71]
    # instances = [25]
    for i in instances
        @testset "Instance: p$(i)" begin
            # Load UFLP data
            data = read_uflp_benchmark_data("p$(i)")
            
            # Solve using standard MIP model
            milp = create_milp(data)
            set_optimizer(milp.model, CPLEX.Optimizer)
            optimize!(milp.model)
            mip_obj = objective_value(milp.model)
            
            # Solve using different Benders decomposition strategies
            loop_strategy = Sequential()
            cut_strategies = Dict(
                "Classical" => ClassicalCut(),
                "FatKnapsack" => FatKnapsackCut()
            )
            params = BendersParams(
                600.0,
                0.00001,
                solver,
                Dict("solver" => solver),
                Dict("solver" => solver),
                Dict(),
                false
            )
            benders_UB = Dict()
            benders_LB = Dict()
            for (name, strategy) in cut_strategies
                result = run_Benders(data, loop_strategy, strategy, params)
                benders_LB[name] = result[end, :LB]
                benders_UB[name] = result[end, :UB]
                @test isapprox(mip_obj, result[end, :LB], atol=0.1)
                @test isapprox(mip_obj, result[end, :UB], atol=0.1)
            end
            @printf("Instance: p%d | %-5s: %.2f | ", i, "MIP", mip_obj)
            for (cut_type, lb) in benders_LB
                @printf("%-5s: %.2f | ", string(cut_type), lb)
            end
            println()
        end
    end
end
