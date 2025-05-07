# Include SCFLP model files
include("$(dirname(dirname(@__DIR__)))/example/scflp/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/scflp/model.jl")

@testset verbose = true "SCFLP Sequential Benders Tests" begin
    # Specify instances to test
    instances = 1:5  # For quick testing
    
    for i in instances
        @testset "Instance: f25-c50-s64-r10-$i" begin
            @info "Testing SCFLP instance $i"
            
            # Load problem data
            problem = read_stochastic_capacited_facility_location_problem("f25-c50-s64-r10-$i")
            
            benders_param = BendersBnBParam(;
                time_limit = 200.0,
                gap_tolerance = 1e-6,
                verbose = false
            )
            
            # Common solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPXPARAM_Threads" => 4)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 4)
            typical_oracle_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)

            # Initialize data object
            dim_x = problem.n_facilities
            dim_t = problem.n_scenarios
            c_x = problem.fixed_costs
            c_t = fill(1/problem.n_scenarios, problem.n_scenarios)
            data = Data(dim_x, dim_t, problem, c_x, c_t)
                
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
                    @info "solving SCFLP f25-c50-s64-r10-$i - classical oracle - no seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)
                    typical_oracle = SeparableOracle(data, ClassicalOracle(), data.problem.n_scenarios; solver_param = typical_oracle_solver_param)
                    for j=1:typical_oracle.N
                        update_model!(typical_oracle.oracles[j], data, j)
                    end
                    root_preprocessing = NoRootNodePreprocessing()
                    lazy_callback = LazyCallback(typical_oracle)
                    user_callback = NoUserCallback()
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
                @testset "Seq" begin
                    @info "solving SCFLP f25-c50-s64-r10-$i - classical oracle - seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)
                    typical_oracle = SeparableOracle(data, ClassicalOracle(), data.problem.n_scenarios; solver_param = typical_oracle_solver_param)
                    for j=1:typical_oracle.N
                        update_model!(typical_oracle.oracles[j], data, j)
                    end
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
                    @info "solving SCFLP f25-c50-s64-r10-$i - classical oracle - seqinout..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)
                    typical_oracle = SeparableOracle(data, ClassicalOracle(), data.problem.n_scenarios; solver_param = typical_oracle_solver_param)
                    for j=1:typical_oracle.N
                        update_model!(typical_oracle.oracles[j], data, j)
                    end
                    root_seq_type = BendersSeqInOut
                    root_param = BendersSeqInOutParam(
                        time_limit = 300.0,
                        gap_tolerance = 1e-6,
                        stabilizing_x = ones(data.dim_x),
                        α = 0.9,
                        λ = 0.1,
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
            
            # Test knapsack oracle
            @testset "Knapsack oracle" begin
                @testset "NoSeq" begin
                    @info "solving SCFLP f25-c50-s64-r10-$i - knapsack oracle - no seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)
                    typical_oracle = SeparableOracle(data, CFLKnapsackOracle(), data.problem.n_scenarios; solver_param = typical_oracle_solver_param)
                    for j=1:typical_oracle.N
                        update_model!(typical_oracle.oracles[j], data, j)
                    end
                    root_preprocessing = NoRootNodePreprocessing()
                    lazy_callback = LazyCallback(typical_oracle)
                    user_callback = NoUserCallback()
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
                @testset "Seq" begin
                    @info "solving SCFLP f25-c50-s64-r10-$i - knapsack oracle - seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)
                    typical_oracle = SeparableOracle(data, CFLKnapsackOracle(), data.problem.n_scenarios; solver_param = typical_oracle_solver_param)
                    for j=1:typical_oracle.N
                        update_model!(typical_oracle.oracles[j], data, j)
                    end
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
                    @info "solving SCFLP f25-c50-s64-r10-$i - knapsack oracle - seqinout..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)
                    typical_oracle = SeparableOracle(data, CFLKnapsackOracle(), data.problem.n_scenarios; solver_param = typical_oracle_solver_param)
                    for j=1:typical_oracle.N
                        update_model!(typical_oracle.oracles[j], data, j)
                    end
                    root_seq_type = BendersSeqInOut
                    root_param = BendersSeqInOutParam(
                        time_limit = 300.0,
                        gap_tolerance = 1e-6,
                        stabilizing_x = ones(data.dim_x),
                        α = 0.9,
                        λ = 0.1,
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
        end
    end
end