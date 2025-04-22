using Test
using JuMP
using Gurobi, CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition
import BendersDecomposition: generate_cuts

include("$(dirname(@__DIR__))/example/uflp/data_reader.jl")
include("$(dirname(@__DIR__))/example/uflp/oracle.jl")
include("$(dirname(@__DIR__))/example/uflp/model.jl")

@testset verbose = true "UFLP Sequential Benders Tests -- MIP master" begin
    # instances = setdiff(1:71, [67])
    instances = [49] # only instance that uses SpecializedBendersSeq
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
            specialized_benders_param = SpecializedBendersSeqParam(;
                                        time_limit = 200.0,
                                        gap_tolerance = 1e-6,
                                        integrality_tolerance = 1e-4,
                                        verbose = true
                                    )
            dcglp_param = DcglpParam(;
                                    time_limit = 1000.0, 
                                    gap_tolerance = 1e-3, 
                                    halt_limit = 250, 
                                    iter_limit = 250,
                                    verbose = true
                            )
            # solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9)
            typical_oracal_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)
            dcglp_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)
                            
            # solve mip for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Default Benders decomposition" begin
                
            end

            @testset "Classical oracle" begin
                @testset "SpecialSeq" begin        
                    @info "solving p$i - classical oracle - seq..."
                    for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf], disjunctive_cut_append_rule in [DisjunctiveCutsSmallerIndices()]
                    # for strengthened in [false], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0], disjunctive_cut_append_rule in [DisjunctiveCutsSmallerIndices()]
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule" begin
                            master = Master(data; solver_param = master_solver_param)
                            update_model!(master, data)

                            typical_oracles = [ClassicalOracle(data; solver_param = typical_oracal_solver_param); ClassicalOracle(data; solver_param = typical_oracal_solver_param)] # for kappa & nu
                            for k=1:2
                                update_model!(typical_oracles[k], data)
                            end

                            # define disjunctive_oracle_attributes and add all the setting to there
                            disjunctive_oracle = DisjunctiveOracle(data, typical_oracles; 
                                                                   solver_param = dcglp_solver_param,
                                                                   param = dcglp_param) 
                            oracle_param = DisjunctiveOracleParam(norm = LpNorm(p), 
                                                                    split_index_selection_rule = LargestFractional(),
                                                                    disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
                                                                    strengthened=strengthened, 
                                                                    add_benders_cuts_to_master=add_benders_cuts_to_master, 
                                                                    fraction_of_benders_cuts_to_master = 0.5, 
                                                                    reuse_dcglp=reuse_dcglp)
                            set_parameter!(disjunctive_oracle, oracle_param)
                            update_model!(disjunctive_oracle, data)

                            env = SpecializedBendersSeq(data, master, disjunctive_oracle; param = specialized_benders_param)
                            
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)

                        end
                    end
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
                # fat-knapsack-based disjunctive cut has sparse gamma_t, so adding only disjunctive cut does not improve lower bound, setting add_benders_cuts_to_master = true
                @testset "Seq" begin
                    @info "solving p$i - fat Knapsack oracle - seq..."
                    for strengthened in [true; false], add_benders_cuts_to_master in [true], reuse_dcglp in [true; false], p in [1.0; Inf], disjunctive_cut_append_rule in [DisjunctiveCutsSmallerIndices()]
                    # for strengthened in [false], add_benders_cuts_to_master in [false], reuse_dcglp in [false], p in [2.0], disjunctive_cut_append_rule in [DisjunctiveCutsSmallerIndices()]
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule" begin
                            master = Master(data; solver_param = master_solver_param)
                            update_model!(master, data)

                            # model-free knapsack-based cuts
                            typical_oracles = [UFLKnapsackOracle(data); UFLKnapsackOracle(data)] # for kappa & nu

                            disjunctive_oracle = DisjunctiveOracle(data, typical_oracles; 
                                                                   solver_param = dcglp_solver_param,
                                                                   param = dcglp_param) 
                            oracle_param = DisjunctiveOracleParam(norm = LpNorm(p), 
                                                                    split_index_selection_rule = LargestFractional(),
                                                                    disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
                                                                    strengthened=strengthened, 
                                                                    add_benders_cuts_to_master=add_benders_cuts_to_master, 
                                                                    fraction_of_benders_cuts_to_master = 0.5, 
                                                                    reuse_dcglp=reuse_dcglp)
                            # norm is used in the initialization.
                            set_parameter!(disjunctive_oracle, oracle_param)
                            update_model!(disjunctive_oracle, data)
                            
                            env = SpecializedBendersSeq(data, master, disjunctive_oracle; param = specialized_benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end
                    end
                end
            end

        end
    end
end