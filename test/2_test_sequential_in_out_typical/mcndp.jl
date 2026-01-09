include("$(dirname(dirname(@__DIR__)))/example/mcndp/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/mcndp/model.jl")


@testset verbose = true "MCNDP Sequential Benders Tests" begin
    # for instance in [0], snipno in [0], budget in [30.0]
        # @testset "instance $instance; snipno $snipno budget $budget" begin
            problem = read_mcndp_instance("R", "", "r01.1.dow")
            # problem = read_mcndp_instance("LargeScaleInstances", "ClassA", "08_05a_.dow")

            # initialize dim_x, dim_t, c_x, c_t
            dim_x = problem.num_arcs
            dim_t = 1
            c_x = problem.fixed_costs
            c_t = [1]
            data = Data(dim_x, dim_t, problem, c_x, c_t)

            # loop parameters
            benders_inout_param = BendersSeqInOutParam(;
            time_limit = 200.0,
            gap_tolerance = 1e-6,
            verbose = true,
            stabilizing_x = ones(data.dim_x),
            α = 0.9,
            λ = 0.1
            )

            # solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPXPARAM_Threads" => 7)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 7, "CPX_PARAM_SCRIND" => 0)
            typical_oracle_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_SCRIND" => 0, "CPX_PARAM_PREIND" => 0)
                        
            # oracle parameters & corepoint
            rtol, atol = 1e-9, 1e-9
            # core_point = fill(data.problem.budget/length(data.problem.D)-1e-3, dim_x)

            # solve mip for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)
            @info mip_opt_val

            @testset "Classic oracle" begin
                # @info "solving CFLP p$i - classical oracle - seq..."
                master = Master(data; solver_param = master_solver_param)
                update_model!(master, data)

                # Construct oracle and set parameters
                classical_param = ClassicalOracleParam(rtol = rtol, atol = atol)
                oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param, oracle_param = classical_param)
                update_model!(oracle, data)

                env = BendersSeqInOut(data, master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end 
        # end
    # end
end