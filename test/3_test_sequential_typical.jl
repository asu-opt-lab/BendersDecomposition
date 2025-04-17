# To-Do: 
# 1. need to be able to change the setting for SeqInOut: e.g., stabilizing point
# 2. assign attributes to the structure, not to JuMP Model: e.g., one may want to setting for CFLKnapsackOracle (e.g., slim, add_only_violated_cuts)
# Done:
# slim cut should be averaged, instead of summation, for numerical stability
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

            # loop parameters
            benders_param = BendersSeqParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = true
                        )
            # solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9)
            typical_oracal_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)

            # solve mip for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Classic oracle" begin
                @testset "SeqInOut" begin
                    @info "solving p$i - classical oracle - seqInOut..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    oracle = ClassicalOracle(data; solver_param = typical_oracal_solver_param)
                    update_model!(oracle, data)
                    
                    stabilizing_x = ones(data.dim_x)
                    env = BendersSeqInOut(data, master, oracle, stabilizing_x; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    # if env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                    # elseif env.termination_status == TimeLimit()
                    #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                    #     @test env.log.LB <= mip_opt_val <= env.log.UB
                    # elseif env.termination_status == InfeasibleOrNumericalIssue()
                    #     @test false
                    # end
                end
                
                @testset "Seq" begin        
                    @info "solving p$i - classical oracle - seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    oracle = ClassicalOracle(data; solver_param = typical_oracal_solver_param)
                    update_model!(oracle, data)

                    env = BendersSeq(data, master, oracle; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    # if env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                    # elseif env.termination_status == TimeLimit()
                    #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                    #     @test env.log.LB <= mip_opt_val <= env.log.UB
                    # elseif env.termination_status == InfeasibleOrNumericalIssue()
                    #     @test false
                    # end
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
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    # model-free knapsack-based cuts
                    oracle = UFLKnapsackOracle(data) 
                    set_parameter!(oracle, "add_only_violated_cuts", true)

                    env = BendersSeq(data, master, oracle; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    # if env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                    # elseif env.termination_status == TimeLimit()
                    #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                    #     @test env.log.LB <= mip_opt_val <= env.log.UB
                    # elseif env.termination_status == InfeasibleOrNumericalIssue()
                    #     @test false
                    # end
                end
                @testset "SeqInOut" begin
                    @info "solving p$i - fat Knapsack oracle - seqInOut..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    # model-free knapsack-based cuts
                    oracle = UFLKnapsackOracle(data) 
                    set_parameter!(oracle, "add_only_violated_cuts", true)

                    stabilizing_x = ones(data.dim_x)
                    env = BendersSeqInOut(data, master, oracle, stabilizing_x; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    # if env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                    # elseif env.termination_status == TimeLimit()
                    #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                    #     @test env.log.LB <= mip_opt_val <= env.log.UB
                    # elseif env.termination_status == InfeasibleOrNumericalIssue()
                    #     @test false
                    # end
                end 
            end

            @testset "slim knapsack oracle" begin
                @testset "Seq" begin
                    @info "solving p$i - slim Knapsack oracle - seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    # model-free knapsack-based cuts
                    oracle = UFLKnapsackOracle(data) # add_only_violated_cuts = true makes it very slow
                    set_parameter!(oracle, "add_only_violated_cuts", false)
                    set_parameter!(oracle, "slim", true)

                    env = BendersSeq(data, master, oracle; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    # if env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                    # elseif env.termination_status == TimeLimit()
                    #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                    #     @test env.log.LB <= mip_opt_val <= env.log.UB
                    # elseif env.termination_status == InfeasibleOrNumericalIssue()
                    #     @test false
                    # end
                end
                @testset "SeqInOut" begin
                    @info "solving p$i - slim Knapsack oracle - seqInOut..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    # model-free knapsack-based cuts
                    oracle = UFLKnapsackOracle(data) 
                    set_parameter!(oracle, "add_only_violated_cuts", false)
                    set_parameter!(oracle, "slim", true)

                    stabilizing_x = ones(data.dim_x)
                    env = BendersSeqInOut(data, master, oracle, stabilizing_x; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    # if env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                    # elseif env.termination_status == TimeLimit()
                    #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                    #     @test env.log.LB <= mip_opt_val <= env.log.UB
                    # elseif env.termination_status == InfeasibleOrNumericalIssue()
                    #     @test false
                    # end
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

            # loop parameters
            benders_param = BendersSeqParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = true
                        )
            # solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9)
            typical_oracal_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)

            # solve mip for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Classic oracle" begin
                @testset "SeqInOut" begin
                    @info "solving p$i - classical oracle - seqInOut..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    oracle = ClassicalOracle(data; solver_param = typical_oracal_solver_param)
                    update_model!(oracle, data)

                    stabilizing_x = ones(data.dim_x)
                    env = BendersSeqInOut(data, master, oracle, stabilizing_x; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    # if env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                    # elseif env.termination_status == TimeLimit()
                    #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                    #     @test env.log.LB <= mip_opt_val <= env.log.UB
                    # elseif env.termination_status == InfeasibleOrNumericalIssue()
                    #     @test false
                    # end
                end
                @testset "Seq" begin        
                    @info "solving p$i - classical oracle - seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    oracle = ClassicalOracle(data; solver_param = typical_oracal_solver_param)
                    update_model!(oracle, data)

                    env = BendersSeq(data, master, oracle; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    # if env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                    # elseif env.termination_status == TimeLimit()
                    #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                    #     @test env.log.LB <= mip_opt_val <= env.log.UB
                    # elseif env.termination_status == InfeasibleOrNumericalIssue()
                    #     @test false
                    # end
                end
            end 
            @testset "Knapsack oracle" begin
                @testset "Seq" begin
                    @info "solving p$i - knapsack oracle - seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    oracle = CFLKnapsackOracle(data; solver_param = typical_oracal_solver_param)
                    update_model!(oracle, data)

                    env = BendersSeq(data, master, oracle; param = benders_param)
                    log = solve!(env)
                    # if env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                    # elseif env.termination_status == TimeLimit()
                    #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                    #     @test env.log.LB <= mip_opt_val <= env.log.UB
                    # elseif env.termination_status == InfeasibleOrNumericalIssue()
                    #     @test false
                    # end
                end
                @testset "SeqInOut" begin
                    @info "solving p$i - knapsack oracle - seqInOut..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    oracle = CFLKnapsackOracle(data; solver_param = typical_oracal_solver_param)
                    update_model!(oracle, data)

                    stabilizing_x = ones(data.dim_x)
                    env = BendersSeqInOut(data, master, oracle, stabilizing_x; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    # if env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                    # elseif env.termination_status == TimeLimit()
                    #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                    #     @test env.log.LB <= mip_opt_val <= env.log.UB
                    # elseif env.termination_status == InfeasibleOrNumericalIssue()
                    #     @test false
                    # end
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

            # loop parameters
            benders_param = BendersSeqParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = true
                        )
            # solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9)
            typical_oracal_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)

            # solve mip for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Classic oracle" begin
                @testset "SeqInOut" begin
                    @info "solving f25-c50-s64-r10-$i - classical oracle - seqInOut..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    oracle = SeparableOracle(data, ClassicalOracle(), data.problem.n_scenarios; solver_param = typical_oracal_solver_param)
                    for j=1:oracle.N
                        update_model!(oracle.oracles[j], data, j)
                    end

                    stabilizing_x = ones(data.dim_x)
                    env = BendersSeqInOut(data, master, oracle, stabilizing_x; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    # if env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                    # elseif env.termination_status == TimeLimit()
                    #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                    #     @test env.log.LB <= mip_opt_val <= env.log.UB
                    # elseif env.termination_status == InfeasibleOrNumericalIssue()
                    #     @test false
                    # end
                end
                
                @testset "Seq" begin        
                    @info "solving f25-c50-s64-r10-$i - classical oracle - seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    oracle = SeparableOracle(data, ClassicalOracle(), data.problem.n_scenarios; solver_param = typical_oracal_solver_param)
                    for j=1:oracle.N
                        update_model!(oracle.oracles[j], data, j)
                    end

                    env = BendersSeq(data, master, oracle; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    # if env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                    # elseif env.termination_status == TimeLimit()
                    #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                    #     @test env.log.LB <= mip_opt_val <= env.log.UB
                    # elseif env.termination_status == InfeasibleOrNumericalIssue()
                    #     @test false
                    # end
                end
            end 
            @testset "Knapsack oracle" begin
                @testset "SeqInOut" begin
                    @info "solving f25-c50-s64-r10-$i - knapsack oracle - seqInOut..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    oracle = SeparableOracle(data, CFLKnapsackOracle(), data.problem.n_scenarios; solver_param = typical_oracal_solver_param)
                    for j=1:oracle.N
                        update_model!(oracle.oracles[j], data, j)
                    end

                    stabilizing_x = ones(data.dim_x)
                    env = BendersSeqInOut(data, master, oracle, stabilizing_x; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    # if env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                    # elseif env.termination_status == TimeLimit()
                    #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                    #     @test env.log.LB <= mip_opt_val <= env.log.UB
                    # elseif env.termination_status == InfeasibleOrNumericalIssue()
                    #     @test false
                    # end
                end
                
                @testset "Seq" begin        
                    @info "solving f25-c50-s64-r10-$i - knapsack oracle - seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    oracle = SeparableOracle(data, CFLKnapsackOracle(), data.problem.n_scenarios; solver_param = typical_oracal_solver_param)
                    for j=1:oracle.N
                        update_model!(oracle.oracles[j], data, j)
                    end

                    env = BendersSeq(data, master, oracle; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    # if env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                    # elseif env.termination_status == TimeLimit()
                    #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                    #     @test env.log.LB <= mip_opt_val <= env.log.UB
                    # elseif env.termination_status == InfeasibleOrNumericalIssue()
                    #     @test false
                    # end
                end
            end
        end
    end
end