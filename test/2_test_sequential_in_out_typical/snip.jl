include("$(dirname(dirname(@__DIR__)))/example/snip/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/snip/model.jl")

@testset verbose = true "SNIP Sequential Benders Tests" begin
    for instance in [0], snipno in [0], budget in [30.0]
        @testset "instance $instance; snipno $snipno budget $budget" begin
            # Load problem data if necessary
            problem = read_snip_data(instance, snipno, budget)

            # initialize dim_x, dim_t, c_x, c_t
            dim_x = length(problem.D)
            dim_t = problem.num_scenarios
            c_x = zeros(dim_x)
            c_t = map(p -> p[3], problem.scenarios)
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
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPXPARAM_Threads" => 4, "CPX_PARAM_SCRIND" => 0)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 4, "CPX_PARAM_SCRIND" => 0)
            typical_oracle_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_SCRIND" => 0)
            basic_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_SCRIND" => 0)

            # oracle parameters & corepoint
            rtol, atol = 1e-9, 1e-9
            core_point = fill(data.problem.budget/length(data.problem.D)-1e-3, dim_x)

            # solve mip for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Classic oracle" begin     
                @info "solving SNIP instance$instance snipno $snipno budget $budget - classical oracle - seqinout..."
                master = Master(data; solver_param = master_solver_param)
                update_model!(master, data)

                # Construct oracle and set parameters
                classical_param = ClassicalOracleParam(rtol = rtol, atol = atol) 
                oracle = SeparableOracle(data, ClassicalOracle(), data.problem.num_scenarios; solver_param = typical_oracle_solver_param, sub_oracle_param = classical_param)
                for j=1:oracle.N
                    update_model!(oracle.oracles[j], data, j)
                end

                env = BendersSeqInOut(data, master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end 

            @testset "Pareto oracle" begin     
                @info "solving SNIP instance$instance snipno $snipno budget $budget - pareto oracle - seqinout..."
                master = Master(data; solver_param = master_solver_param)
                update_model!(master, data)
                
                # Construct oracle and set parameters
                pareto_param = ParetoOracleParam(rtol = rtol, atol = atol, core_point = core_point)
                oracle = SeparableOracle(data, ParetoOracle(), data.problem.num_scenarios; solver_param = basic_solver_param, sub_oracle_param = pareto_param)
        
                for j=1:oracle.N
                    update_model!(oracle.oracles[j], data, j)
                    model_reformulation!(oracle.oracles[j])
                end

                env = BendersSeqInOut(data, master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end 

            @testset "Unified oracle" begin
                @info "solving SNIP instance$instance snipno $snipno budget $budget - unified oracle - seqinout..."
                master = Master(data; solver_param = master_solver_param)
                update_model!(master, data)

                unified_param = UnifiedOracleParam(rtol = rtol, atol = atol)
                oracle = SeparableOracle(data, UnifiedOracle(), data.problem.num_scenarios; solver_param = typical_oracle_solver_param, sub_oracle_param = unified_param)
                for j=1:oracle.N
                    update_model!(oracle.oracles[j], data, j)
                    model_reformulation!(oracle.oracles[j])
                end

                env = BendersSeqInOut(data, master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end
        end
    end
end