# to be overwritten, they should be included outside testset
include("$(dirname(dirname(@__DIR__)))/example/scflp/data_reader.jl")
# include("$(dirname(@__DIR__))/example/cflp/oracle.jl")
include("$(dirname(dirname(@__DIR__)))/example/scflp/model.jl")

@testset verbose = true "Stochastic CFLP Sequential In/Out Benders Tests" begin
    # instances = setdiff(1:71, [67])
    instances = 1:5
    for i in instances
        @testset "Instance: f25-c50-s64-r10-$i" begin
            # Load problem data if necessary
            problem = read_stochastic_capacited_facility_location_problem("f25-c50-s64-r10-$i")
            
            # initialize dim_x, dim_t, c_x, c_t
            dim_x = problem.n_facilities
            dim_t = problem.n_scenarios
            c_x = problem.fixed_costs
            c_t = fill(1/problem.n_scenarios, problem.n_scenarios)
            data = Data(dim_x, dim_t, problem, c_x, c_t)
            @assert dim_x == length(data.c_x)
            @assert dim_t == length(data.c_t)

            # loop parameters
            benders_inout_param = BendersSeqInOutParam(;
                time_limit = 200.0,
                gap_tolerance = 1e-6,
                verbose = false,
                stabilizing_x = ones(data.dim_x),
                α = 0.9,
                λ = 0.1
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
                @info "solving SCFLP f25-c50-s64-r10-$i - classical oracle - seqInOut..."
                master = Master(data; solver_param = master_solver_param)
                update_model!(master, data)

                oracle = SeparableOracle(data, ClassicalOracle(), data.problem.n_scenarios; solver_param = typical_oracal_solver_param)
                for j=1:oracle.N
                    update_model!(oracle.oracles[j], data, j)
                end

                env = BendersSeqInOut(data, master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end 
            
            @testset "Knapsack oracle" begin
                @info "solving SCFLP f25-c50-s64-r10-$i - knapsack oracle - seqInOut..."
                master = Master(data; solver_param = master_solver_param)
                update_model!(master, data)

                oracle = SeparableOracle(data, CFLKnapsackOracle(), data.problem.n_scenarios; solver_param = typical_oracal_solver_param)
                for j=1:oracle.N
                    update_model!(oracle.oracles[j], data, j)
                end

                env = BendersSeqInOut(data, master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end
        end
    end
end