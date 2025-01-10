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
    # uflp_files = readdir("src/Datasets/locssall", join=true)
    solver = "CPLEX"
    for i in [1:66; 68:71]
        @testset "Instance: p$(i)" begin
            # Load UFLP data
            data = read_uflp_benchmark_data("p$(i)")
            
            # Solve using standard MIP model
            milp = create_milp(data)
            set_optimizer(milp.model, CPLEX.Optimizer)
            optimize!(milp.model)
            mip_objective = objective_value(milp.model)
            
            # Solve using different Benders decomposition strategies
            loop_strategy = Sequential()
            cut_strategies = [ClassicalCut(), FatKnapsackCut(), SlimKnapsackCut()]
            benders_objectives = Dict()
            
            for cut_strategy in cut_strategies
                master = create_master_problem(data, cut_strategy)
                set_optimizer(master.model, CPLEX.Optimizer)
                sub = create_sub_problem(data, cut_strategy)
                if cut_strategy == ClassicalCut()
                    set_optimizer(sub.model, CPLEX.Optimizer)
                end
                params = BendersParams(7200.0, 1e-9, solver, Dict("solver" => solver), Dict("solver" => solver), Dict("solver" => solver), true) 

                result = run_Benders(data, loop_strategy, cut_strategy, params)
                benders_objectives[typeof(cut_strategy)] = result[end, :UB]

                # Compare results
                @test isapprox(mip_objective, benders_objectives[typeof(cut_strategy)], rtol=1e-4)
            end
            
            # Print results for all strategies
            @printf("Instance: p%d | %-5s: %.4f | ", i, "MIP", mip_objective)
            for (cut_type, objective) in benders_objectives
                @printf("%-5s: %.4f | ", string(cut_type), objective)
            end
            println()
        end
    end
end
