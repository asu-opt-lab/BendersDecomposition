using BendersDecomposition
using Test
using JuMP

@testset verbose = true "Stochastic CFLP Sequential In/Out Benders Tests" begin
    instances = 1:5

    for i in instances
        @testset "Instance: f25-c50-s64-r10-$i" begin
            # Load problem data
            problem = read_stochastic_capacited_facility_location_problem("f25-c50-s64-r10-$i")
            
            # Initialize data object
            dim_x = problem.n_facilities
            dim_t = problem.n_scenarios
            c_x = problem.fixed_costs
            c_t = fill(1/problem.n_scenarios, problem.n_scenarios)
            
            data = Data(dim_x, dim_t, problem, c_x, c_t)
            @assert dim_x == length(data.c_x)
            @assert dim_t == length(data.c_t)

            # Loop parameters
            benders_inout_param = BendersSeqInOutParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = false,
                            stabilizing_x = ones(data.dim_x),
                            α = 0.9,
                            λ = 0.1
                        )
            
            # Solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6)
            master_solver_param = Dict("solver" => "CPLEX", "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6)
            typical_oracle_solver_param = Dict("solver" => "CPLEX", "CPXPARAM_Threads" => 7, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1)

            # Solve MIP for reference
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

                oracle = SeparableOracle(data, ClassicalOracle(), data.problem.n_scenarios; solver_param = typical_oracle_solver_param)
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

                oracle = SeparableOracle(data, CFLKnapsackOracle(), data.problem.n_scenarios; solver_param = typical_oracle_solver_param)
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