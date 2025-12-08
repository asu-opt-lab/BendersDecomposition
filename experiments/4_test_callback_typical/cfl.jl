using BendersDecomposition
using Test
using JuMP

@testset verbose = true "CFLP Callback Benders Tests" begin
    instances = setdiff(1:71, [67])
    
    for i in instances
        @testset "Instance: p$i" begin
            @info "Testing CFLP easy instance $i"
            
            # Load problem data
            problem = read_cflp_benchmark_data("p$i")

            # BnB parameters
            benders_param = BendersBnBParam(;
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
                @testset "NoSeq" begin
                    @info "solving CFLP p$i - classical oracle - no seq..."
                    master = Master(problem; customize = customize_master_model!)
                    oracle = ClassicalOracle(problem, master; customize = customize_sub_model!)

                    root_preprocessing = NoRootNodePreprocessing()
                    lazy_callback = LazyCallback(oracle)
                    user_callback = NoUserCallback()

                    env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end

                @testset "Seq" begin
                    @info "solving CFLP p$i - classical oracle - seq..."
                    master = Master(problem; customize = customize_master_model!)
                    oracle = ClassicalOracle(problem, master; customize = customize_sub_model!)

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
                    @info "solving CFLP p$i - classical oracle - seqinout..."
                    master = Master(problem; customize = customize_master_model!)
                    oracle = ClassicalOracle(problem, master; customize = customize_sub_model!)

                    root_seq_type = BendersSeqInOut
                    root_param = BendersSeqInOutParam(
                                time_limit = 300.0,
                                gap_tolerance = 1e-9,
                                stabilizing_x = ones(problem.n_facilities),
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
            
            @testset "Knapsack oracle" begin
                @testset "NoSeq" begin
                    @info "solving CFLP p$i - knapsack oracle - no seq..."
                    master = Master(problem; customize = customize_master_model!)
                    oracle = CFLKnapsackOracle(problem, master; customize = customize_sub_model!)

                    root_preprocessing = NoRootNodePreprocessing()
                    lazy_callback = LazyCallback(oracle)
                    user_callback = NoUserCallback()

                    env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
                
                @testset "Seq" begin
                    @info "solving CFLP p$i - knapsack oracle - seq..."
                    master = Master(problem; customize = customize_master_model!)
                    oracle = CFLKnapsackOracle(problem, master; customize = customize_sub_model!)

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
                    @info "solving CFLP p$i - knapsack oracle - seqinout..."
                    master = Master(problem; customize = customize_master_model!)
                    oracle = CFLKnapsackOracle(problem, master; customize = customize_sub_model!)

                    root_seq_type = BendersSeqInOut
                    root_param = BendersSeqInOutParam(
                                time_limit = 300.0,
                                gap_tolerance = 1e-9,
                                stabilizing_x = ones(problem.n_facilities),
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

    # instances = 1:5
    # for i in instances
    #     @testset "Instance: p$i" begin
    #         @info "Testing CFLP hard instance $i"
            
    #         # Load problem data
    #         problem = read_GK_data("f100-c100-r5-$i")
            
    #         # Get standard parameters
    #         benders_param, dcglp_param, mip_solver_param, master_solver_param, 
    #         typical_oracle_solver_param, dcglp_solver_param = get_standard_params()
            
    #         # Create data object
    #         data = create_data(problem)
            
    #         # Solve MIP for reference
    #         mip_opt_val = solve_reference_mip(data, mip_solver_param)
            
    #         # Standard test info
    #         test_info = (mip_opt_val, master_solver_param, typical_oracle_solver_param)
            
    #         # Test classical oracle
    #         @testset "Classic oracle" begin
    #             oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param)
    #             update_model!(oracle, data)
    #             run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "classical oracle")
    #         end
            
    #         # Test CFLKnapsack oracle
    #         @testset "CFLKnapsack oracle" begin
    #             oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
    #             update_model!(oracle, data)
    #             run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "CFLKnapsack oracle")
    #         end
    #     end
    # end
end