using BendersX
using Test
using JuMP
using CPLEX

@testset verbose = true "UFLP Callback Benders Tests" begin
    instances = setdiff(1:71, [67])
    
    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data
            data = read_uflp_benchmark_data("p$i")

            # BnB parameters
            benders_param = BendersBnBParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = false
                        )

            # Solve MIP for reference
            mip_model = Model()
            customize_mip_model!(mip_model, data)
            set_optimizer_attribute(mip_model, "CPX_PARAM_BRDIR", 1)
            optimize!(mip_model)
            @assert termination_status(mip_model) == OPTIMAL
            mip_opt_val = objective_value(mip_model)
            
            @testset "Classic oracle" begin
                @testset "NoSeq" begin
                    @info "solving UFLP p$i - classical oracle - no seq..."
                    # This setting can use default initializer
                    master = Master(data; customize = customize_master_model!)
                    set_optimizer_attribute(master.model, "CPX_PARAM_BRDIR", 1)
                    oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
                
                    # root_preprocessing = NoRootNodePreprocessing()
                    # lazy_callback = LazyCallback(oracle)
                    # user_callback = NoUserCallback()
                    # env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)

                    env = BendersBnB(master, oracle; param = benders_param)

                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end

                @testset "Seq" begin
                    @info "solving UFLP p$i - classical oracle - seq..."
                    
                    master = Master(data; customize = customize_master_model!)
                    set_optimizer_attribute(master.model, "CPX_PARAM_BRDIR", 1)
                    oracle = ClassicalOracle(data, master; customize = customize_sub_model!)

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
                    @info "solving UFLP p$i - classical oracle - seqinout..."
                    master = Master(data; customize = customize_master_model!)
                    set_optimizer_attribute(master.model, "CPX_PARAM_BRDIR", 1)
                    oracle = ClassicalOracle(data, master; customize = customize_sub_model!)

                    root_seq_type = BendersSeqInOut
                    root_param = BendersSeqInOutParam(
                                time_limit = 300.0,
                                gap_tolerance = 1e-9,
                                stabilizing_x = ones(data.n_facilities),
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

                function customize_master_model!(model::Model, data::UFLPData)
                    optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_BRDIR" => 1, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)
                    set_optimizer(model, optimizer)
                    @variable(model, x[1:data.n_facilities], Bin)
                    @variable(model, t[1:data.n_customers] >= -1e6)
                    @constraint(model, sum(x) >= 2)
                    @objective(model, Min, data.fixed_costs'* x + sum(t))
                    return (x = x, ), t
                end

                @testset "NoSeq" begin
                    @info "solving UFLP p$i - fat knapsack oracle - no seq..."
                    # This setting can use default initializer
                    master = Master(data; customize = customize_master_model!)
                    oracle = UFLKnapsackOracle(data) 
                    set_parameter!(oracle, "add_only_violated_cuts", true)
                    
                    # root_preprocessing = NoRootNodePreprocessing()
                    # lazy_callback = LazyCallback(oracle)
                    # user_callback = NoUserCallback()
                    # env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)

                    env = BendersBnB(master, oracle; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end

                @testset "Seq" begin
                    @info "solving UFLP p$i - fat knapsack oracle - seq..."
                    master = Master(data; customize = customize_master_model!)
                    oracle = UFLKnapsackOracle(data) 
                    set_parameter!(oracle, "add_only_violated_cuts", true)

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
                    @info "solving UFLP p$i - fat knapsack oracle - seqinout..."
                    master = Master(data; customize = customize_master_model!)
                    oracle = UFLKnapsackOracle(data) 
                    set_parameter!(oracle, "add_only_violated_cuts", true)

                    root_seq_type = BendersSeqInOut
                    root_param = BendersSeqInOutParam(
                                time_limit = 300.0,
                                gap_tolerance = 1e-9,
                                stabilizing_x = ones(data.n_facilities),
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
            # To test slim version, users can use # set_parameter!(oracle, "slim", true)
        end
    end
end