include("$(dirname(dirname(@__DIR__)))/example/uflp/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/uflp/oracle.jl")
include("$(dirname(dirname(@__DIR__)))/example/uflp/model.jl")

@testset verbose = true "UFLP Sequential Benders Tests" begin
    instances = setdiff(1:71, [67])
    # instances = 30:35
    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data if necessary
            problem = read_uflp_benchmark_data("p$(i)")
            
            # initialize dim_x, dim_t, c_x, c_t
            dim_x = problem.n_facilities
            c_x = problem.fixed_costs
            dim_t = 1 # classical cut
            c_t = [1]
            
            data = Data(dim_x, dim_t, problem, c_x, c_t)
            @assert dim_x == length(data.c_x)
            @assert dim_t == length(data.c_t)

            # loop parameters
            benders_param = BendersSeqParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = false
                        )
            # solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPXPARAM_Threads" => 4)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 4)
            typical_oracal_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)

            # solve mip for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Classic oracle" begin
                @info "solving UFLP p$i - classical oracle - seq..."
                master = Master(data; solver_param = master_solver_param)
                update_model!(master, data)

                oracle = ClassicalOracle(data; solver_param = typical_oracal_solver_param)
                update_model!(oracle, data)

                env = BendersSeq(data, master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end

            # initialize dim_x, dim_t, c_x, c_t
            dim_x = problem.n_facilities
            c_x = problem.fixed_costs
            dim_t = problem.n_customers # knapsack cut
            c_t = ones(dim_t)
            
            data = Data(dim_x, dim_t, problem, c_x, c_t)
            @assert dim_x == length(data.c_x)
            @assert dim_t == length(data.c_t)

            @testset "fat knapsack oracle" begin
                @info "solving UFLP p$i - fat Knapsack oracle - seq..."
                master = Master(data; solver_param = master_solver_param)
                update_model!(master, data)

                # model-free knapsack-based cuts
                oracle = UFLKnapsackOracle(data) 
                set_parameter!(oracle, "add_only_violated_cuts", true)

                env = BendersSeq(data, master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end

            # @testset "slim knapsack oracle" begin
            #     @info "solving UFLP p$i - slim Knapsack oracle - seq..."
            #     master = Master(data; solver_param = master_solver_param)
            #     update_model!(master, data)

            #     # model-free knapsack-based cuts
            #     oracle = UFLKnapsackOracle(data) # add_only_violated_cuts = true makes it very slow
            #     set_parameter!(oracle, "add_only_violated_cuts", false)
            #     set_parameter!(oracle, "slim", true)

            #     env = BendersSeq(data, master, oracle; param = benders_param)
            #     log = solve!(env)
            #     @test env.termination_status == Optimal()
            #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            # end
        end
    end
end