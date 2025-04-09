# To-Do:
# when no fractional value exists, randomly choose any index
# streamline the definition for attributes: model vs dcglp
# verbose for master and oracle directly
# add_benders_cuts_to_master only for violated ones by (x_value, t_value)
# test strengthening technique -- should work
# test multiple scenarios -- should work
# test fat knapsack -- should work
# add disjunctive cut handler for dcglp
# add lifting
# even if disjunctiveOracle does not generate a cut, we should not terminate it; dependent on the split 
# define BendersSeq and BendersSeqInOut and BendersBnB
# define BendersSeq and BendersSeqInOut and BendersBnB for disjunctive? continuous master?

# Issues:
# 1. DisjunctiveOracle: Setting CPX_PARAM_EPRHS tightly (e.g., < 1e-6) results in dcglp terminating with ALMOST_INFEASIBLE --> set it tightly (1e-9) and outputs a typical Benders cut when dcglp is ALMOST_INFEASIBLE
# 2. DisjunctiveOracle: Setting zero_tol in solve_dcglp! large (e.g., 1e-4) results in disjunctive Benders with reuse_dcglp = true terminating with incorrect solution --> set it tightly as 1e-9
# 3. DisjunctiveOracle: solve_dcglp! becomes stall since the true violation should be multipled with omega_value[:z], which can be fairly small --> terminate dcglp when LB does not improve for a certain number of iterations
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

@testset verbose = true "UFLP Sequential Benders Tests" begin
    instances = setdiff(1:71, [67])
    # instances = 21:30
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
                    disjunctive_oracle = DisjunctiveOracle(data, typical_oracles, LpNorm(Inf), RandomFractional(); strengthened=true, add_benders_cuts_to_master=false, reuse_dcglp=true, verbose=true) # when not reusing, all combinations are correct
                    assign_attributes!(disjunctive_oracle.dcglp, params.oracle_attributes)

                    env = BendersEnv(data, master, disjunctive_oracle, Seq())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end
                # # SeqInOut for disjunctive oracle is actually unnecessary extra, but it works though
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
                #     disjunctive_oracle = DisjunctiveOracle(data, typical_oracles, LpNorm(1.0), RandomFractional(); strengthened=true, add_benders_cuts_to_master=false, reuse_dcglp=true, verbose=true)
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

            # @testset "fat knapsack oracle" begin
            #     @testset "Seq" begin
            #         @info "solving p$i - fat Knapsack oracle - seq..."
            #         master = Master(data)
            #         assign_attributes!(master.model, params.master_attributes)
            #         # problem specific constraints
            #         update_model!(master, data)

            #         # model-free knapsack-based cuts
            #         typical_oracles = [UFLKnapsackOracle(data); UFLKnapsackOracle(data)] # for kappa & nu

            #         # define disjunctive_oracle_attributes and add all the setting to there
            #         disjunctive_oracle = DisjunctiveOracle(data, typical_oracles, LpNorm(1.0), RandomFractional(); strengthened=false, add_benders_cuts_to_master=false, reuse_dcglp=false, verbose=true)
            #         # streamline the definition for attributes: model vs dcglp
            #         assign_attributes!(disjunctive_oracle.dcglp, params.oracle_attributes)

            #         env = BendersEnv(data, master, disjunctive_oracle, Seq())
            #         run_Benders(env, params)
            #         @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
            #     end
            # end
            @testset "slim knapsack oracle" begin
                @testset "Seq" begin
                    @info "solving p$i - slim Knapsack oracle - seq..."
                    master = Master(data)
                    assign_attributes!(master.model, params.master_attributes)
                    # problem specific constraints
                    update_model!(master, data)

                    # model-free knapsack-based cuts
                    typical_oracles = [UFLKnapsackOracle(data, slim=true); UFLKnapsackOracle(data, slim=true)] # for kappa & nu

                    # define disjunctive_oracle_attributes and add all the setting to there
                    disjunctive_oracle = DisjunctiveOracle(data, typical_oracles, LpNorm(1.0), RandomFractional(); strengthened=true, add_benders_cuts_to_master=false, reuse_dcglp=false, verbose=true)
                    # streamline the definition for attributes: model vs dcglp
                    assign_attributes!(disjunctive_oracle.dcglp, params.oracle_attributes)

                    env = BendersEnv(data, master, disjunctive_oracle, Seq())
                    run_Benders(env, params)
                    @test isapprox(mip_opt_val, env.master.obj_value, atol=1e-5)
                end
            end
        end
    end
end