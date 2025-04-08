using Test
using JuMP
using CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition
import BendersDecomposition: generate_cuts
# global_logger(ConsoleLogger(stderr, Logging.Warn))

include("$(dirname(@__DIR__))/example/uflp/data_reader.jl")
include("$(dirname(@__DIR__))/example/uflp/oracle.jl")
include("$(dirname(@__DIR__))/example/uflp/model.jl")

@testset verbose = true "UFLP Sequential Benders Tests" begin
    # instances = setdiff(1:71, [67])
    instances = 30:35
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

            params = BendersParams(
                200.0,
                0.00001,
                Dict("solver" => "CPLEX"),
                Dict("solver" => "CPLEX"),
                true
            )

            # solve mip for reference
            mip = Mip(data)
            assign_attributes!(mip.model, params.master_attributes)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Classic oracle" begin
                @testset "SeqInOut" begin
                    @info "solving p$i - classical oracle - seqInOut..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    oracle = ClassicalOracle(data)
                    assign_attributes!(oracle.model, params.oracle_attributes)
                    # problem specific 
                    update_model!(oracle, data)

                    env = BendersEnv(data, master, oracle, SeqInOut())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end
                
                @testset "Seq" begin        
                    @info "solving p$i - classical oracle - seq..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    oracle = ClassicalOracle(data)
                    assign_attributes!(oracle.model, params.oracle_attributes)
                    # problem specific 
                    update_model!(oracle, data)

                    env = BendersEnv(data, master, oracle, Seq())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end
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
                @testset "Seq" begin
                    @info "solving p$i - fat Knapsack oracle - seq..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    # model-free knapsack-based cuts
                    oracle = UFLKnapsackOracle(data) 

                    env = BendersEnv(data, master, oracle, Seq())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end
                @testset "SeqInOut" begin
                    @info "solving p$i - fat Knapsack oracle - seqInOut..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    # model-free knapsack-based cuts
                    oracle = UFLKnapsackOracle(data) 

                    env = BendersEnv(data, master, oracle, SeqInOut())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end 
            end

            @testset "slim knapsack oracle" begin
                ## interestingly, aggregating only violated cuts slows down Benders significantly
                @testset "Seq" begin
                    @info "solving p$i - slim Knapsack oracle - seq..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    # model-free knapsack-based cuts
                    oracle = UFLKnapsackOracle(data; slim=true) 

                    env = BendersEnv(data, master, oracle, Seq())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end
                @testset "SeqInOut" begin
                    @info "solving p$i - slim Knapsack oracle - seqInOut..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    # model-free knapsack-based cuts
                    oracle = UFLKnapsackOracle(data; slim=true) 

                    env = BendersEnv(data, master, oracle, SeqInOut())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end 
            end
        end
    end
end

# to be overwritten, they should be included outside testset
include("$(dirname(@__DIR__))/example/cflp/data_reader.jl")
include("$(dirname(@__DIR__))/example/cflp/oracle.jl")
include("$(dirname(@__DIR__))/example/cflp/model.jl")

@testset verbose = true "CFLP Sequential Benders Tests" begin
    # instances = setdiff(1:71, [67])
    instances = 30:35
    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data if necessary
            problem = read_cflp_benchmark_data("p$i")
            # problem = read_GK_data("f100-c100-r3-1")
            
            # initialize dim_x, dim_t, c_x, c_t
            dim_x = problem.n_facilities
            dim_t = 1
            c_x = problem.fixed_costs
            c_t = [1]
            data = Data(dim_x, dim_t, problem, c_x, c_t)
            @assert dim_x == length(data.c_x)
            @assert dim_t == length(data.c_t)

            params = BendersParams(
                200.0,
                0.00001,
                Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9),
                Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9),
                false
            )

            # solve mip for reference
            mip = Mip(data)
            assign_attributes!(mip.model, params.master_attributes)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Classic oracle" begin
                @testset "SeqInOut" begin
                    @info "solving p$i - classical oracle - seqInOut..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    oracle = ClassicalOracle(data)
                    assign_attributes!(oracle.model, params.oracle_attributes)
                    # problem specific 
                    update_model!(oracle, data)

                    env = BendersEnv(data, master, oracle, SeqInOut())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end
                
                @testset "Seq" begin        
                    @info "solving p$i - classical oracle - seq..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    oracle = ClassicalOracle(data)
                    assign_attributes!(oracle.model, params.oracle_attributes)
                    # problem specific 
                    update_model!(oracle, data)

                    env = BendersEnv(data, master, oracle, Seq())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end
            end 
            @testset "Knapsack oracle" begin
                @testset "Seq" begin
                    @info "solving p$i - knapsack oracle - seq..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    oracle = CFLKnapsackOracle(data)
                    assign_attributes!(oracle.model, params.oracle_attributes)
                    # problem specific 
                    update_model!(oracle, data)

                    env = BendersEnv(data, master, oracle, Seq())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end
                @testset "SeqInOut" begin
                    @info "solving p$i - knapsack oracle - seqInOut..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    oracle = CFLKnapsackOracle(data)
                    assign_attributes!(oracle.model, params.oracle_attributes)
                    # problem specific 
                    update_model!(oracle, data)

                    env = BendersEnv(data, master, oracle, SeqInOut())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end 
            end
        end
    end
end

# to be overwritten, they should be included outside testset
include("$(dirname(@__DIR__))/example/scflp/data_reader.jl")
# include("$(dirname(@__DIR__))/example/cflp/oracle.jl")
include("$(dirname(@__DIR__))/example/scflp/model.jl")

@testset verbose = true "Stochastic CFLP Sequential Benders Tests" begin
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

            params = BendersParams(
                200.0,
                0.00001,
                Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9),
                Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9),
                false
            )

            # solve mip for reference
            mip = Mip(data)
            assign_attributes!(mip.model, params.master_attributes)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Classic oracle" begin
                @testset "SeqInOut" begin
                    @info "solving f25-c50-s64-r10-$i - classical oracle - seqInOut..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    oracle = SeparableOracle(data, ClassicalOracle(), data.problem.n_scenarios)
                    for j=1:oracle.N
                        assign_attributes!(oracle.oracles[j].model, params.oracle_attributes)
                        # problem specific 
                        update_model!(oracle.oracles[j], data, j)
                    end

                    env = BendersEnv(data, master, oracle, SeqInOut())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end
                
                @testset "Seq" begin        
                    @info "solving f25-c50-s64-r10-$i - classical oracle - seq..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    oracle = SeparableOracle(data, ClassicalOracle(), data.problem.n_scenarios)
                    for j=1:oracle.N
                        assign_attributes!(oracle.oracles[j].model, params.oracle_attributes)
                        # problem specific 
                        update_model!(oracle.oracles[j], data, j)
                    end

                    env = BendersEnv(data, master, oracle, Seq())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end
            end 
            @testset "Knapsack oracle" begin
                @testset "SeqInOut" begin
                    @info "solving f25-c50-s64-r10-$i - classical oracle - seqInOut..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    oracle = SeparableOracle(data, CFLKnapsackOracle(), data.problem.n_scenarios)
                    for j=1:oracle.N
                        assign_attributes!(oracle.oracles[j].model, params.oracle_attributes)
                        # problem specific 
                        update_model!(oracle.oracles[j], data, j)
                    end

                    env = BendersEnv(data, master, oracle, SeqInOut())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end
                
                @testset "Seq" begin        
                    @info "solving f25-c50-s64-r10-$i - classical oracle - seq..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    oracle = SeparableOracle(data, CFLKnapsackOracle(), data.problem.n_scenarios)
                    for j=1:oracle.N
                        assign_attributes!(oracle.oracles[j].model, params.oracle_attributes)
                        # problem specific 
                        update_model!(oracle.oracles[j], data, j)
                    end

                    env = BendersEnv(data, master, oracle, Seq())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end
            end
        end
    end
end