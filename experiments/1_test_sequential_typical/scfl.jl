using BendersDecomposition
using Test
using JuMP

@testset verbose = true "Stochastic CFLP Sequential Benders Tests" begin
    instances = 1:5
    
    for i in instances
        @testset "Instance: f25-c50-s64-r10-$i" begin
            
            # Load problem data
            problem = read_stochastic_capacited_facility_location_problem("f25-c50-s64-r10-$i")
            
            # Loop parameters
            benders_param = BendersSeqParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = false
                        )
            
            # Solve MIP for reference
            mip_model = Model()
            customize_mip_model!(mip_model, problem)
            optimize!(mip_model)
            @assert termination_status(mip_model) == OPTIMAL
            mip_opt_val = objective_value(mip_model)            
            
            @testset "Classic oracle" begin     
                @info "solving SCFLP f25-c50-s64-r10-$i - classical oracle - seq..."
                master = Master(problem; customize = customize_master_model!)
                oracle = SeparableOracle(problem, master, ClassicalOracle(), problem.n_scenarios; customize = customize_sub_model!)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end 
            
            @testset "Knapsack oracle" begin
                @info "solving SCFLP f25-c50-s64-r10-$i - knapsack oracle - seq..."
                master = Master(problem; customize = customize_master_model!)
                oracle = SeparableOracle(problem, master, CFLKnapsackOracle(), problem.n_scenarios; customize = customize_sub_model!)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end
        end
    end
end