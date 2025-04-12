# To-Do:
# add X to dcglp
# streamline the definition for attributes: (1) model vs dcglp; (2) termination parameter for dcglp; (3) verbose for master and oracle directly
# test multiple scenarios -- should work
# add lifting
# continuous master; check split index rule for continuous master
# define BendersSeq and BendersSeqInOut and BendersBnB
# define BendersSeq and BendersSeqInOut and BendersBnB for disjunctive? continuous master?
# plug in other solver, e.g., conic solvers like Mosek 

# Done:
# when no fractional value exists, randomly choose any index --> done; need to check select_disjunctive_inequality for LP master
# even if disjunctiveOracle does not generate a cut, we should not terminate it --> done; when latest tau is close to zero, return typical Benders cut
# test fat knapsack --> works
# test strengthening technique --> works
# add_benders_cuts_to_master only for a farction of violated ones by (x_value, t_value) --> done
# add disjunctive cut handler for dcglp

# Issues:
# 1. DisjunctiveOracle: Setting CPX_PARAM_EPRHS tightly (e.g., < 1e-6) results in dcglp terminating with ALMOST_INFEASIBLE --> set it tightly (1e-9) and outputs a typical Benders cut when dcglp is ALMOST_INFEASIBLE
# 2. DisjunctiveOracle: Setting zero_tol in solve_dcglp! large (e.g., 1e-4) results in disjunctive Benders with reuse_dcglp = true terminating with incorrect solution --> set it tightly as 1e-9
# 3. DisjunctiveOracle: solve_dcglp! becomes stall since the true violation should be multipled with omega_value[:z], which can be fairly small --> terminate dcglp when LB does not improve for a certain number of iterations
# 4. DisjunctiveOracle: the fat-knapsack-based disjunctive cut may have a sparse gamma_t, so adding only disjunctive cut does not improve lower bound, add_benders_cuts_to_master should be set at true
using Test
using JuMP
using Gurobi, CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition
import BendersDecomposition: generate_cuts
# global_logger(ConsoleLogger(stderr, Logging.Warn))

include("$(dirname(@__DIR__))/example/uflp/data_reader.jl")
include("$(dirname(@__DIR__))/example/uflp/oracle.jl")
include("$(dirname(@__DIR__))/example/uflp/model.jl")

@testset verbose = true "UFLP Sequential Benders Tests -- MIP master" begin
    # instances = setdiff(1:71, [67])
    instances = 21:21
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
                200.0, # time_limit
                1e-6,
                Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9),
                Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9), 
                true
            )

            # solve mip for reference
            mip = Mip(data)
            assign_attributes!(mip.model, params.master_attributes)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Classical oracle" begin
                @testset "Seq" begin        
                    @info "solving p$i - classical oracle - seq..."
                    for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf], disjunctive_cut_append_rule in [NoDisjunctiveCuts(); AllDisjunctiveCuts(); DisjunctiveCutsSmallerIndices()]
                        # for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [Inf], disjunctive_cut_append_rule in [DisjunctiveCutsSmallerIndices()]
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule" begin
                            master = Master(data)
                            assign_attributes!(master.model, params.master_attributes)
                            # problem specific constraints
                            update_model!(master, data)

                            typical_oracles = [ClassicalOracle(data); ClassicalOracle(data)] # for kappa & nu
                            for k=1:2
                                assign_attributes!(typical_oracles[k].model, params.oracle_attributes)
                                # problem specific 
                                update_model!(typical_oracles[k], data)
                            end

                            # define disjunctive_oracle_attributes and add all the setting to there
                            disjunctive_oracle = DisjunctiveOracle(data, typical_oracles, LpNorm(p), RandomFractional(), disjunctive_cut_append_rule; strengthened=strengthened, add_benders_cuts_to_master=add_benders_cuts_to_master, fraction_of_benders_cuts_to_master = 0.5, reuse_dcglp=reuse_dcglp, verbose=false) # when not reusing, all combinations are correct
                            update_model!(disjunctive_oracle, data)
                            assign_attributes!(disjunctive_oracle.dcglp, params.oracle_attributes)

                            env = BendersEnv(data, master, disjunctive_oracle, Seq())
                            run_Benders(env, params)
                            @test env.log.termination_status == Optimal()
                            # if env.log.termination_status == Optimal()
                                @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                            # elseif env.log.termination_status == TimeLimit()
                            #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                            #     @test env.log.LB <= mip_opt_val <= env.log.UB
                            # elseif env.log.termination_status == InfeasibleOrNumericalIssue()
                            #     @test false
                            # end
                        end
                    end
                end
                # # SeqInOut for DisjunctiveOracle is actually unnecessary extra, but it works though
                # @testset "SeqInOut" begin
                #     @info "solving p$i - classical oracle - seqInOut..."
                #     master = Master(data)
                #     assign_attributes!(master.model, params.master_attributes)
                #     # problem specific constraints
                #     update_model!(master, data)

                #     typical_oracles = [ClassicalOracle(data); ClassicalOracle(data)] # for kappa & nu
                #     for k=1:2
                #         assign_attributes!(typical_oracles[k].model, params.oracle_attributes)
                #         # problem specific 
                #         update_model!(typical_oracles[k], data)
                #     end
                #     # define disjunctive_oracle_attributes and add all the setting to there
                #     disjunctive_oracle = DisjunctiveOracle(data, typical_oracles, LpNorm(1.0), RandomFractional(); strengthened=true, add_benders_cuts_to_master=false, reuse_dcglp=true, verbose=false)
                #     # streamline the definition for attributes: model vs dcglp
                #     assign_attributes!(disjunctive_oracle.dcglp, params.oracle_attributes)

                #     env = BendersEnv(data, master, disjunctive_oracle, SeqInOut())
                #     run_Benders(env, params)
                #     @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                # end
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
                    for strengthened in [true; false], add_benders_cuts_to_master in [true], reuse_dcglp in [true; false], p in [1.0; Inf], disjunctive_cut_append_rule in [NoDisjunctiveCuts(); AllDisjunctiveCuts(); DisjunctiveCutsSmallerIndices()]
                        # for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [Inf], disjunctive_cut_append_rule in [DisjunctiveCutsSmallerIndices()]
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule" begin
                            master = Master(data)
                            assign_attributes!(master.model, params.master_attributes)
                            # problem specific constraints
                            update_model!(master, data)

                            # model-free knapsack-based cuts
                            typical_oracles = [UFLKnapsackOracle(data); UFLKnapsackOracle(data)] # for kappa & nu

                            # define disjunctive_oracle_attributes and add all the setting to there
                            disjunctive_oracle = DisjunctiveOracle(data, typical_oracles, LpNorm(p), RandomFractional(), disjunctive_cut_append_rule; strengthened=strengthened, add_benders_cuts_to_master=add_benders_cuts_to_master, fraction_of_benders_cuts_to_master = 0.5, reuse_dcglp=reuse_dcglp, verbose=false)
                            update_model!(disjunctive_oracle, data)
                            assign_attributes!(disjunctive_oracle.dcglp, params.oracle_attributes)

                            env = BendersEnv(data, master, disjunctive_oracle, Seq())
                            run_Benders(env, params)
                            @test env.log.termination_status == Optimal()
                            # if env.log.termination_status == Optimal()
                                @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                            # elseif env.log.termination_status == TimeLimit()
                            #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                            #     @test env.log.LB <= mip_opt_val <= env.log.UB
                            # elseif env.log.termination_status == InfeasibleOrNumericalIssue()
                            #     @test false
                            # end
                        end
                    end
                end
            end
            @testset "slim knapsack oracle" begin
                @testset "Seq" begin
                    @info "solving p$i - slim Knapsack oracle - seq..."
                    for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf], disjunctive_cut_append_rule in [NoDisjunctiveCuts(); AllDisjunctiveCuts(); DisjunctiveCutsSmallerIndices()]
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule" begin
                            master = Master(data)
                            assign_attributes!(master.model, params.master_attributes)
                            # problem specific constraints
                            update_model!(master, data)

                            # model-free knapsack-based cuts
                            typical_oracles = [UFLKnapsackOracle(data, slim=true); UFLKnapsackOracle(data, slim=true)] # for kappa & nu

                            # define disjunctive_oracle_attributes and add all the setting to there
                            disjunctive_oracle = DisjunctiveOracle(data, typical_oracles, LpNorm(1.0), RandomFractional(), AllDisjunctiveCuts(); strengthened=strengthened, add_benders_cuts_to_master=add_benders_cuts_to_master, fraction_of_benders_cuts_to_master = 0.5, reuse_dcglp=reuse_dcglp, verbose=false)
                            # streamline the definition for attributes: model vs dcglp
                            update_model!(disjunctive_oracle, data)
                            assign_attributes!(disjunctive_oracle.dcglp, params.oracle_attributes)

                            env = BendersEnv(data, master, disjunctive_oracle, Seq())
                            run_Benders(env, params)
                            @test env.log.termination_status == Optimal()
                            # if env.log.termination_status == Optimal()
                                @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                            # elseif env.log.termination_status == TimeLimit()
                            #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                            #     @test env.log.LB <= mip_opt_val <= env.log.UB
                            # elseif env.log.termination_status == InfeasibleOrNumericalIssue()
                            #     @test false
                            # end
                        end
                    end
                end
            end
        end
    end
end