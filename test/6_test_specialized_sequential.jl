# tight tolerance for dcglp; adjust tol for generate_cut with omega_0
# dcglp is not solved -> cannot proceed the algorithm, and should be terminated. Added an optimal argument to dcglp `throw_typical_cuts_for_errors=false`
# remove other_constraints from classical oracle.
# add a method for generating optimal vertex for SpecializedBendersSeq

using Test
using JuMP
using Gurobi, CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition
import BendersDecomposition: generate_cuts
global_logger(ConsoleLogger(stderr, Logging.Debug))

# loop parameters
specialized_benders_param = SpecializedBendersSeqParam(;
                                                        time_limit = 1000.0,
                                                        lp_gap_tolerance = 1e-9,
                                                        integrality_tolerance = 1e-9,
                                                        verbose = true
                                                        )
dcglp_param = DcglpParam(;
                        time_limit = 1000.0, 
                        gap_tolerance = 1e-9, 
                        halt_limit = Int(1e9), 
                        iter_limit = Int(1e9),
                        verbose = true
                        )
# solver parameters
mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9)
# master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_LPMETHOD" => 1, "CPX_PARAM_EPOPT" => 1e-9)
master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_LPMETHOD" => 1, "CPX_PARAM_EPOPT" => 1e-9)
typical_oracal_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)
dcglp_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_LPMETHOD" => 3)

include("$(dirname(@__DIR__))/example/uflp/data_reader.jl")
include("$(dirname(@__DIR__))/example/uflp/oracle.jl")
include("$(dirname(@__DIR__))/example/uflp/model.jl")

@testset verbose = true "UFLP Specialized Sequential Benders Tests" begin
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
                            
            # solve mip for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            @testset "Classical oracle" begin
                @testset "SpecialSeq" begin        
                    # for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf]
                    # for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [Inf] # too slow

                    for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0]
                    # for strengthened in [false], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0]
                    # for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [1.0]
                    # for strengthened in [false], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [1.0]

                    # for strengthened in [false], add_benders_cuts_to_master in [false], reuse_dcglp in [false], p in [1.0]
                    # for strengthened in [true], add_benders_cuts_to_master in [false], reuse_dcglp in [false], p in [1.0]
                    # for strengthened in [true], add_benders_cuts_to_master in [false], reuse_dcglp in [true], p in [1.0] # fail
                    # for strengthened in [false], add_benders_cuts_to_master in [false], reuse_dcglp in [true], p in [1.0] #fail
                        @info "solving p$i - begin oracle - Specialized seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p"
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p" begin
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
                                                                    disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(), 
                                                                    strengthened=strengthened, 
                                                                    add_benders_cuts_to_master=add_benders_cuts_to_master, 
                                                                    fraction_of_benders_cuts_to_master = 0.5, 
                                                                    reuse_dcglp=reuse_dcglp)
                            set_parameter!(disjunctive_oracle, oracle_param)
                            update_model!(disjunctive_oracle, data)

                            env = SpecializedBendersSeq(data, master, disjunctive_oracle; param = specialized_benders_param)
                            
                            log = solve!(env)
                            @test env.termination_status == Optimal() ? isapprox(mip_opt_val, env.obj_value, atol=1e-5) : false

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
                @testset "SpecialSeq" begin
                    for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0] 
                    # for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf] 
                        @info "solving p$i - fat Knapsack oracle - Specialized seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p"
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p" begin
                            master = Master(data; solver_param = master_solver_param)
                            update_model!(master, data)

                            # model-free knapsack-based cuts
                            typical_oracles = [UFLKnapsackOracle(data); UFLKnapsackOracle(data)] # for kappa & nu

                            disjunctive_oracle = DisjunctiveOracle(data, typical_oracles; 
                                                                   solver_param = dcglp_solver_param,
                                                                   param = dcglp_param) 
                            oracle_param = DisjunctiveOracleParam(norm = LpNorm(p), 
                                                                    split_index_selection_rule = LargestFractional(),
                                                                    disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(), 
                                                                    strengthened=strengthened, 
                                                                    add_benders_cuts_to_master=add_benders_cuts_to_master, 
                                                                    fraction_of_benders_cuts_to_master = 0.5, 
                                                                    reuse_dcglp=reuse_dcglp)
                            # norm is used in the initialization.
                            set_parameter!(disjunctive_oracle, oracle_param)
                            update_model!(disjunctive_oracle, data)
                            
                            env = SpecializedBendersSeq(data, master, disjunctive_oracle; param = specialized_benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal() ? isapprox(mip_opt_val, env.obj_value, atol=1e-5) : false
                        end
                    end
                end
            end
        end
    end
end

# to be overwritten, they should be included outside testset
include("$(dirname(@__DIR__))/example/cflp/data_reader.jl")
include("$(dirname(@__DIR__))/example/cflp/oracle.jl")
include("$(dirname(@__DIR__))/example/cflp/model.jl")

@testset verbose = true "CFLP Specialized Sequential Benders Tests" begin
    # instances = setdiff(1:71, [67])
    instances = [25 32 34 36 49 51]
    # numerical issue: 29 30 31 35
    # success: 25 32 34 36 49 51
    # lp: 1:24, 26:28, 33, 37-48, 50, 52-66, 68-71
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
                            
            # solve mip for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)
            x_opt = value.(mip.model[:x])
            t_opt = value.(mip.model[:t])
            @info "MIP_opt_val: $mip_opt_val"
            @info "MIP sol: $(x_opt)"

            @testset "Classical oracle" begin
                @testset "SpecialSeq" begin        
                    # for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf]
                    for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0]
                        @info "solving CFLP p$i - disjunctive oracle/classical - specialized seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p"
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p" begin
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
                                                                    disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(), 
                                                                    strengthened=strengthened, 
                                                                    add_benders_cuts_to_master=add_benders_cuts_to_master, 
                                                                    fraction_of_benders_cuts_to_master = 0.5, 
                                                                    reuse_dcglp=reuse_dcglp)
                            set_parameter!(disjunctive_oracle, oracle_param)
                            update_model!(disjunctive_oracle, data)

                            env = SpecializedBendersSeq(data, master, disjunctive_oracle; param = specialized_benders_param)
                            
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            # if env.log.termination_status == Optimal()
                            if !isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                                t_opt_ = [sum(t_opt)]
                                @error "Failed ***** mip_opt_val = $(mip_opt_val) vs BD_obj_val = $(env.obj_value)" 
                                # optimize!(master.model)
                                opt_sol = Dict{VariableRef, Float64}()
                                for i = 1:data.dim_x
                                    opt_sol[master.model[:x][i]] = x_opt[i]
                                end
                                for i = 1:data.dim_t
                                    opt_sol[master.model[:t][i]] = t_opt_[i]
                                end
                                
                                @info primal_feasibility_report(env.master.model, opt_sol)
                                @info data.c_x' * x_opt + data.c_t' * t_opt_
                                
                                for v in keys(opt_sol)
                                    fix(v, opt_sol[v]; force=true)
                                end
                                optimize!(env.master.model)
                                @info objective_value(env.master.model)
                            end
                            @test env.termination_status == Optimal() ? isapprox(mip_opt_val, env.obj_value, atol=1e-5) : false
                        end
                    end
                end
            end 

            @testset "Knapsack oracle" begin
                @testset "SpecialSeq" begin
                    # for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf] 
                    for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0] 
                        @info "solving CFLP p$i - disjunctive oracle/knapsack- specialized seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p"
                            @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p" begin
                            master = Master(data; solver_param = master_solver_param)
                            update_model!(master, data)

                            typical_oracles = [CFLKnapsackOracle(data; solver_param = typical_oracal_solver_param); CFLKnapsackOracle(data; solver_param = typical_oracal_solver_param)]
                            for k=1:2
                                update_model!(typical_oracles[k], data)
                            end

                            disjunctive_oracle = DisjunctiveOracle(data, typical_oracles; 
                                                                   solver_param = dcglp_solver_param,
                                                                   param = dcglp_param) 
                            oracle_param = DisjunctiveOracleParam(norm = LpNorm(p), 
                                                                    split_index_selection_rule = LargestFractional(),
                                                                    disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(), 
                                                                    strengthened=strengthened, 
                                                                    add_benders_cuts_to_master=add_benders_cuts_to_master, 
                                                                    fraction_of_benders_cuts_to_master = 0.5, 
                                                                    reuse_dcglp=reuse_dcglp)
                            # norm is used in the initialization.
                            set_parameter!(disjunctive_oracle, oracle_param)
                            update_model!(disjunctive_oracle, data)
                            
                            env = SpecializedBendersSeq(data, master, disjunctive_oracle; param = specialized_benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal() ? isapprox(mip_opt_val, env.obj_value, atol=1e-5) : false
                        end
                    end
                end
            end
        end
    end
end