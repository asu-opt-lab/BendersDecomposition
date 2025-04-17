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

# to be overwritten, they should be included outside testset
include("$(dirname(@__DIR__))/example/cflp/data_reader.jl")
include("$(dirname(@__DIR__))/example/cflp/oracle.jl")
include("$(dirname(@__DIR__))/example/cflp/model.jl")

@testset verbose = true "CFLP Callback Benders Tests" begin
    # instances = setdiff(1:71, [67])
    instances = 1:1
    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data if necessary
            # problem = read_cflp_benchmark_data("p$i")
            problem = read_GK_data("f100-c100-r5-$i")
            
            # initialize dim_x, dim_t, c_x, c_t
            dim_x = problem.n_facilities
            dim_t = 1
            c_x = problem.fixed_costs
            c_t = [1]
            data = Data(dim_x, dim_t, problem, c_x, c_t)
            @assert dim_x == length(data.c_x)
            @assert dim_t == length(data.c_t)

            # loop parameters
            benders_param = BendersBnBParam(;
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
                # @testset "no root preprocessing" begin        
                #     @info "solving p$i - classical oracle - seq..."
                #     master = Master(data; solver_param = master_solver_param)
                #     update_model!(master, data)

                #     oracle = ClassicalOracle(data; solver_param = typical_oracal_solver_param)
                #     update_model!(oracle, data)

                #     callback_param = BendersCallbackParam(;
                #             time_limit = 400.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true,
                #             preprocessing_type = nothing,
                #             root_param = nothing
                #         )
                #     env = BendersCallback(data, master, oracle; param = callback_param)
                #     log = solve!(env)
                #     @test env.termination_status == Optimal()
                #     # if env.termination_status == Optimal()
                #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                #     # elseif env.termination_status == TimeLimit()
                #     #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                #     #     @test env.log.LB <= mip_opt_val <= env.log.UB
                #     # elseif env.termination_status == InfeasibleOrNumericalIssue()
                #     #     @test false
                #     # end
                # end
                @testset "Seq root preprocessing" begin        
                    @info "solving p$i - classical oracle - seq..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    oracle = ClassicalOracle(data; solver_param = typical_oracal_solver_param)
                    update_model!(oracle, data)

                    root_param = BendersSeqParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = true
                        )
                    callback_param = BendersBnBParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = true,
                            preprocessing_type = BendersSeq,
                            root_param = root_param
                        )
                    lazy_callback = LazyCallbackParam(default_lazy_callback, Dict{String, Any}())
                    user_callback = UserCallbackParam(default_user_callback, Dict{String, Any}("frequency" => 250))
                    env = BendersBnB(data, master, oracle, lazy_callback, user_callback; param = callback_param)
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
                # @testset "SeqInOut" begin
                #     @info "solving p$i - classical oracle - seqInOut..."
                #     master = Master(data; solver_param = master_solver_param)
                #     update_model!(master, data)

                #     oracle = ClassicalOracle(data; solver_param = typical_oracal_solver_param)
                #     update_model!(oracle, data)

                #     stabilizing_x = ones(data.dim_x)
                #     root_param = BendersSeqInOutParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true,
                #             stabilizing_x = stabilizing_x,
                #             α = 0.9,
                #             λ = 0.1
                #         )
                #     callback_param = BendersCallbackParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true,
                #             preprocessing_type = BendersSeqInOut,
                #             root_param = root_param
                #         )
                #     env = BendersCallback(data, master, oracle; param = callback_param)
                #     log = solve!(env)
                #     @test env.termination_status == Optimal()
                #     # if env.termination_status == Optimal()
                #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                #     # elseif env.termination_status == TimeLimit()
                #     #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                #     #     @test env.log.LB <= mip_opt_val <= env.log.UB
                #     # elseif env.termination_status == InfeasibleOrNumericalIssue()
                #     #     @test false
                #     # end
                # end
            end 

            @testset "Knapsack oracle" begin
                # @testset "no root preprocessing" begin        
                #     @info "solving p$i - classical oracle - seq..."
                #     master = Master(data; solver_param = master_solver_param)
                #     update_model!(master, data)

                #     oracle = CFLKnapsackOracle(data; solver_param = typical_oracal_solver_param)
                #     update_model!(oracle, data)

                #     callback_param = BendersCallbackParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true,
                #             preprocessing_type = nothing,
                #             root_param = nothing
                #         )
                #     env = BendersCallback(data, master, oracle; param = callback_param)
                #     log = solve!(env)
                #     @test env.termination_status == Optimal()
                #     # if env.termination_status == Optimal()
                #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                #     # elseif env.termination_status == TimeLimit()
                #     #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                #     #     @test env.log.LB <= mip_opt_val <= env.log.UB
                #     # elseif env.termination_status == InfeasibleOrNumericalIssue()
                #     #     @test false
                #     # end
                # end
                # @testset "Seq root preprocessing" begin        
                #     @info "solving p$i - classical oracle - seq..."
                #     master = Master(data; solver_param = master_solver_param)
                #     update_model!(master, data)

                #     oracle = CFLKnapsackOracle(data; solver_param = typical_oracal_solver_param)
                #     update_model!(oracle, data)

                #     root_param = BendersSeqParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true,
                #         )
                #     callback_param = BendersCallbackParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true,
                #             preprocessing_type = BendersSeq,
                #             root_param = root_param
                #         )
                #     env = BendersCallback(data, master, oracle; param = callback_param)
                #     log = solve!(env)
                #     @test env.termination_status == Optimal()
                #     # if env.termination_status == Optimal()
                #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                #     # elseif env.termination_status == TimeLimit()
                #     #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                #     #     @test env.log.LB <= mip_opt_val <= env.log.UB
                #     # elseif env.termination_status == InfeasibleOrNumericalIssue()
                #     #     @test false
                #     # end
                # end
                # @testset "SeqInOut" begin
                #     @info "solving p$i - classical oracle - seqInOut..."
                #     master = Master(data; solver_param = master_solver_param)
                #     update_model!(master, data)

                #     oracle = CFLKnapsackOracle(data; solver_param = typical_oracal_solver_param)
                #     update_model!(oracle, data)

                #     stabilizing_x = ones(data.dim_x)
                #     root_param = BendersSeqInOutParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true,
                #             stabilizing_x = stabilizing_x,
                #             α = 0.9,
                #             λ = 0.1
                #         )
                #     callback_param = BendersCallbackParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true,
                #             preprocessing_type = BendersSeqInOut,
                #             root_param = root_param
                #         )
                #     env = BendersCallback(data, master, oracle; param = callback_param)
                #     log = solve!(env)
                #     @test env.termination_status == Optimal()
                #     # if env.termination_status == Optimal()
                #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                #     # elseif env.termination_status == TimeLimit()
                #     #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                #     #     @test env.log.LB <= mip_opt_val <= env.log.UB
                #     # elseif env.termination_status == InfeasibleOrNumericalIssue()
                #     #     @test false
                #     # end
                # end
            end 
            
        end
    end
end

