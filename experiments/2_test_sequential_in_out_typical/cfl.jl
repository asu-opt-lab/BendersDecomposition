using BendersDecomposition
using Test
using JuMP

@testset verbose = true "CFLP Sequential In/Out Benders Tests" begin
    instances = setdiff(1:71, [67])

    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data
            data = read_cflp_benchmark_data("p$i")
            
            # Loop parameters
            benders_inout_param = BendersSeqInOutParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = false,
                            stabilizing_x = ones(data.n_facilities),
                            α = 0.9,
                            λ = 0.1
                        )
            
                        # Solve MIP for reference
            mip_model = Model()
            customize_mip_model!(mip_model, data)
            optimize!(mip_model)
            @assert termination_status(mip_model) == OPTIMAL
            mip_opt_val = objective_value(mip_model)

            @testset "Classic oracle" begin
                @info "solving CFLP p$i - classical oracle - seqInOut..."
                master = Master(data; customize = customize_master_model!)
                oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
                env = BendersSeqInOut(master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end 
            
            @testset "Knapsack oracle" begin
                @info "solving CFLP p$i - knapsack oracle - seqInOut..."
                master = Master(data; customize = customize_master_model!)
                oracle = CFLKnapsackOracle(data, master; customize = customize_sub_model!)
                env = BendersSeqInOut(master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end
        end
    end
end