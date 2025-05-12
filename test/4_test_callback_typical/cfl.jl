include("$(dirname(dirname(@__DIR__)))/example/cflp/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/cflp/oracle.jl")
include("$(dirname(dirname(@__DIR__)))/example/cflp/model.jl")

@testset verbose = true "CFLP Callback Benders Tests" begin
    # Specify instances to test
    instances = setdiff(1:71, [67])  # For quick testing
    
    for i in instances
        @testset "Instance: p$i" begin
            @info "Testing CFLP easy instance $i"
            
            # Load problem data
            problem = read_cflp_benchmark_data("p$i")
            
            # Get standard parameters
            benders_param = BendersBnBParam(;
                time_limit = 200.0,
                gap_tolerance = 1e-6,
                verbose = false
            )
            
            # Common solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPXPARAM_Threads" => 4)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 4)
            typical_oracle_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)
            
            # Create data object
            dim_x = problem.n_facilities
            dim_t = 1
            c_x = problem.fixed_costs
            c_t = [1]
            data = Data(dim_x, dim_t, problem, c_x, c_t)
            @assert dim_x == length(data.c_x)
            @assert dim_t == length(data.c_t)
            
            # Solve MIP for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)
            
            # Test classical oracle
            @testset "Classic oracle" begin
                @testset "NoSeq" begin
                    @info "solving CFLP p$i - classical oracle - no seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)
                    typical_oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param)
                    update_model!(typical_oracle, data)
                    root_preprocessing = NoRootNodePreprocessing()
                    lazy_callback = LazyCallback(typical_oracle)
                    user_callback = NoUserCallback()
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
                @testset "Seq" begin
                    @info "solving CFLP p$i - classical oracle - seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)
                    typical_oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param)
                    update_model!(typical_oracle, data)
                    root_seq_type = BendersSeq
                    root_param = BendersSeqParam(;
                        time_limit = 200.0,
                        gap_tolerance = 1e-6,
                        verbose = false
                    )
                    root_preprocessing = RootNodePreprocessing(typical_oracle, root_seq_type, root_param)
                    lazy_callback = LazyCallback(typical_oracle)
                    user_callback = NoUserCallback()
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
                @testset "SeqInOut" begin
                    @info "solving CFLP p$i - classical oracle - seqinout..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)
                    typical_oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param)
                    update_model!(typical_oracle, data)
                    root_seq_type = BendersSeqInOut
                    root_param = BendersSeqInOutParam(
                        time_limit = 300.0,
                        gap_tolerance = 1e-6,
                        stabilizing_x = ones(data.dim_x),
                        α = 0.9,
                        λ = 0.1,
                        verbose = true
                        verbose = false
                    )
                    root_preprocessing = RootNodePreprocessing(typical_oracle, root_seq_type, root_param)
                    lazy_callback = LazyCallback(typical_oracle)
                    user_callback = NoUserCallback()
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
            end
            
            # Test CFLKnapsack oracle
            @testset "CFLKnapsack oracle" begin
                @testset "NoSeq" begin
                    @info "solving CFLP p$i - CFLKnapsack oracle - no seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)
                    oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
                    update_model!(oracle, data)
                    root_preprocessing = NoRootNodePreprocessing()
                    lazy_callback = LazyCallback(oracle)
                    user_callback = NoUserCallback()
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
                @testset "Seq" begin
                    @info "solving CFLP p$i - CFLKnapsack oracle - seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)
                    oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
                    update_model!(oracle, data)
                    root_seq_type = BendersSeq
                    root_param = BendersSeqParam(;
                        time_limit = 200.0,
                        gap_tolerance = 1e-6,
                        verbose = false
                    )
                    root_preprocessing = RootNodePreprocessing(oracle, root_seq_type, root_param)   
                    lazy_callback = LazyCallback(oracle)    
                    user_callback = NoUserCallback()
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end 
                @testset "SeqInOut" begin
                    @info "solving CFLP p$i - CFLKnapsack oracle - seqinout..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)
                    oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
                    update_model!(oracle, data)
                    root_seq_type = BendersSeqInOut
                    root_param = BendersSeqInOutParam(
                        time_limit = 300.0,
                        gap_tolerance = 1e-6,
                        stabilizing_x = ones(data.dim_x),
                        α = 0.9,    
                        λ = 0.1,
                        verbose = false
                    )
                    root_preprocessing = RootNodePreprocessing(oracle, root_seq_type, root_param)
                    lazy_callback = LazyCallback(oracle)
                    user_callback = NoUserCallback()
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
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