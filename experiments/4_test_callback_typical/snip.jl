using BendersX
using Test
using JuMP

@testset verbose = true "SNIP Sequential Benders Tests" begin
    for instance in [0], snipno in [0], budget in [30.0]
        @testset "instance $instance; snipno $snipno budget $budget" begin
            data = read_snip_data(instance, snipno, budget)
            
            # BnB parameters
            benders_param = BendersBnBParam(;
                            time_limit = 200.0,
                            verbose = false
                        )
                
            # Solve MIP for reference
            mip_model = Model()
            customize_mip_model!(mip_model, data)
            optimize!(mip_model)
            @assert termination_status(mip_model) == OPTIMAL
            mip_opt_val = objective_value(mip_model)

            @testset "Classic oracle" begin
                @testset "NoSeq" begin
                    @info "solving SNIP instance$instance snipno $snipno budget $budget - classical oracle - no seq..."
                    master = Master(data; customize = customize_master_model!)
                    oracle = SeparableOracle(data, master, ClassicalOracle(), data.num_scenarios; customize = customize_sub_model!)

                    root_preprocessing = NoRootNodePreprocessing()
                    lazy_callback = LazyCallback(oracle)
                    user_callback = NoUserCallback()

                    env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end

                @testset "Seq" begin
                    @info "solving SNIP instance$instance snipno $snipno budget $budget - classical oracle - seq..."
                    master = Master(data; customize = customize_master_model!)
                    oracle = SeparableOracle(data, master, ClassicalOracle(), data.num_scenarios; customize = customize_sub_model!)

                    root_seq_type = BendersSeq
                    root_param = BendersSeqParam(;
                                time_limit = 200.0,
                                gap_tolerance = 1e-9,
                                verbose = false
                            )

                    root_preprocessing = RootNodePreprocessing(oracle, root_seq_type, root_param)
                    lazy_callback = LazyCallback(oracle)
                    user_callback = NoUserCallback()

                    env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end

                @testset "SeqInOut" begin
                    @info "solving SNIP instance$instance snipno $snipno budget $budget - classical oracle - seqinout..."
                    master = Master(data; customize = customize_master_model!)
                    oracle = SeparableOracle(data, master, ClassicalOracle(), data.num_scenarios; customize = customize_sub_model!)

                    root_seq_type = BendersSeqInOut
                    root_param = BendersSeqInOutParam(
                                time_limit = 300.0,
                                gap_tolerance = 1e-9,
                                stabilizing_x = ones(length(data.D)),
                                α = 0.9,
                                λ = 0.1,
                                verbose = false
                            )

                    root_preprocessing = RootNodePreprocessing(oracle, root_seq_type, root_param)
                    lazy_callback = LazyCallback(oracle)
                    user_callback = NoUserCallback()
                    
                    env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
            end 
        end
    end
end