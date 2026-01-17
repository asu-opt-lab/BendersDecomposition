include("$(dirname(dirname(@__DIR__)))/example/mcndp/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/mcndp/model.jl")

@testset verbose = true "CFLP Callback Benders Tests" begin
    # for i in instances
        # @testset "Instance: p$i" begin
    
            # problem = read_mcndp_instance("C", "", "c35.dow")
            # problem = read_mcndp_instance("R", "", "r11.7.dow")
            problem = read_mcndp_instance("R", "", "r01.7.dow")

            # problem = read_mcndp_instance("LargeScaleInstances", "ClassA", "08_05a_.dow")

            # loop parameters
            benders_param = BendersBnBParam(;
                time_limit = 200.0,
                verbose = true
            )
            
            # Common solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPXPARAM_Threads" => 7, "CPX_PARAM_SCRIND" => 0)

            # Automatic presolve
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_SCRIND" => 0)
            typical_oracle_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_SCRIND" => 0)

            # master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_SCRIND" => 0, "CPX_PARAM_PREIND" => 0)
            # typical_oracle_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_SCRIND" => 0, "CPX_PARAM_PREIND" => 0)

            # master_solver_param = Dict("solver" => "Gurobi", "Presolve" => 2, "OutputFlag" => 0)
            # typical_oracle_solver_param = Dict("solver" => "Gurobi", "Presolve" => 2, "OutputFlag" => 0)

            # master_solver_param = Dict("solver" => "Gurobi", "OutputFlag" => 0)
            # typical_oracle_solver_param = Dict("solver" => "Gurobi", "OutputFlag" => 0)

            # Create data object
            dim_x = problem.num_arcs
            dim_t = 1
            c_x = problem.fixed_costs
            c_t = [1]
            data = Data(dim_x, dim_t, problem, c_x, c_t)

            # oracle parameters & corepoint
            rtol, atol = 1e-6, 1e-6
            
            # Solve MIP for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)
            @info mip_opt_val
            # mip_opt_val = 371475.0
            # mip_opt_val = 2.294912e6

            @testset "Classic oracle" begin
                # @testset "Seq" begin
                #     master = Master(data; solver_param = master_solver_param)
                #     update_model!(master, data)
                #     # Construct oracle and set parameters
                #     classical_param = ClassicalOracleParam(rtol = rtol, atol = atol) 
                #     typical_oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param, oracle_param = classical_param)
                #     update_model!(typical_oracle, data)
                #     root_seq_type = BendersSeq
                #     root_param = BendersSeqParam(;
                #         time_limit = 200.0,
                #         gap_tolerance = 1e-6,
                #         verbose = true
                #     )
                #     root_preprocessing = RootNodePreprocessing(typical_oracle, root_seq_type, root_param)
                #     lazy_callback = LazyCallback(typical_oracle)
                #     user_callback = NoUserCallback()
                #     env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                #     log = solve!(env)
                #     @test env.termination_status == Optimal()
                #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                # end
                @testset "SeqInOut" begin
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)
                    # Construct oracle and set parameters
                    classical_param = ClassicalOracleParam(rtol = rtol, atol = atol) 
                    typical_oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param, oracle_param = classical_param)
                    update_model!(typical_oracle, data)
                    root_seq_type = BendersSeqInOut
                    root_param = BendersSeqInOutParam(
                        time_limit = 300.0,
                        gap_tolerance = 1e-6,
                        stabilizing_x = ones(data.dim_x),
                        α = 0.9,
                        λ = 0.1,
                        verbose = true
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
        # end
    # end
end