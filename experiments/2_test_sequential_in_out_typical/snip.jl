using BendersX
using Test
using JuMP

@testset verbose = true "SNIP Sequential Benders Tests" begin
    for instance in [0], snipno in [0], budget in [30.0]
        @testset "instance $instance; snipno $snipno budget $budget" begin
            # Load problem data
            data = read_snip_data(instance, snipno, budget)

            # Loop parameters
            benders_inout_param = BendersSeqInOutParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = false,
                            stabilizing_x = ones(length(data.D)),
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
                @info "solving SNIP instance-$instance snipno-$snipno budget-$budget - classical oracle - seq..."
                master = Master(data; customize = customize_master_model!)
                oracle = SeparableOracle(data, master, ClassicalOracle(), data.num_scenarios; customize = customize_sub_model!)
                env = BendersSeqInOut(master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end 
        end
    end
end