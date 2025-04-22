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
    instances = 4:4
    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data if necessary
            # problem = read_cflp_benchmark_data("p$i")
            problem = read_GK_data("f700-c700-r5-$i")
            
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
            dcglp_param = DcglpParam(;
                        time_limit = 1000.0, 
                        gap_tolerance = 1e-3, 
                        halt_limit = 3, 
                        iter_limit = 250,
                        verbose = true
                )
            # solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9)
            typical_oracle_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)
            dcglp_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)

            # solve mip for reference
            # mip = Mip(data)
            # assign_attributes!(mip.model, mip_solver_param)
            # update_model!(mip, data)
            # optimize!(mip.model)
            # @assert termination_status(mip.model) == OPTIMAL
            # mip_opt_val = objective_value(mip.model)

            @testset "Classic oracle" begin
                # @testset "No root preprocessing" begin        
                #     @info "solving p$i - classical oracle - no root preprocessing..."
                #     master = Master(data; solver_param = master_solver_param)
                #     update_model!(master, data)

                #     oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param)
                #     update_model!(oracle, data)

                #     root_seq_type = BendersSeq
                #     root_param = BendersSeqParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true
                #         )
                    
                #     # Create root node preprocessing with oracle
                #     root_preprocessing = NoRootNodePreprocessing()
                    
                #     # Create callbacks
                #     lazy_callback = LazyCallback(oracle=oracle)
                #     user_callback = NoUserCallback()
                    
                #     # Create BnB parameter
                #     callback_param = BendersBnBParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true
                #         )
                    
                #     # Create BendersBnB environment
                #     env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param=callback_param)
                    
                #     # Solve
                #     obj_value, elapsed_time = solve!(env)
                    
                #     @test env.termination_status == Optimal()
                #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                # end
                # @testset "Seq root preprocessing" begin        
                #     @info "solving p$i - classical oracle - seq..."
                #     master = Master(data; solver_param = master_solver_param)
                #     update_model!(master, data)

                #     oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param)
                #     update_model!(oracle, data)

                #     root_seq_type = BendersSeq
                #     root_param = BendersSeqParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true
                #         )
                    
                #     # Create root node preprocessing with oracle
                #     root_preprocessing = RootNodePreprocessing(oracle, root_seq_type, root_param)
                    
                #     # Create callbacks
                #     lazy_callback = LazyCallback(oracle=oracle)
                #     user_callback = UserCallback(params=UserCallbackParam(frequency=250), oracle=oracle)
                    
                #     # Create BnB parameter
                #     callback_param = BendersBnBParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true
                #         )
                    
                #     # Create BendersBnB environment
                #     env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param=callback_param)
                    
                #     # Solve
                #     obj_value, elapsed_time = solve!(env)
                    
                #     @test env.termination_status == Optimal()
                #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                # end
                @testset "SeqInOut root preprocessing" begin        
                    @info "solving p$i - classical oracle - seqinout..."
                    master = Master(data; solver_param = master_solver_param)
                    update_model!(master, data)

                    oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param)
                    update_model!(oracle, data)

                    root_seq_type = BendersSeqInOut
                    root_param = BendersSeqInOutParam(;
                            time_limit = 50.0,
                            gap_tolerance = 1e-6,
                            stabilizing_x = ones(data.dim_x),
                            α = 0.9,
                            λ = 0.1,
                            verbose = true
                        )
                    
                    # Create root node preprocessing with oracle
                    root_preprocessing = RootNodePreprocessing(oracle, root_seq_type, root_param)
                    
                    # Create callbacks
                    lazy_callback = LazyCallback(oracle=oracle)
                    user_callback = UserCallback(params=UserCallbackParam(frequency=250), oracle=oracle)
                    
                    # Create BnB parameter
                    callback_param = BendersBnBParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = true
                        )
                    
                    # Create BendersBnB environment
                    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param=callback_param)
                    
                    # Solve
                    obj_value, elapsed_time = solve!(env)
                    
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end
            end

            @testset "CFLKnapsack oracle" begin
                # @testset "No root preprocessing" begin
                #     @info "solving p$i - CFLKnapsack oracle - no root preprocessing..."
                #     master = Master(data; solver_param = master_solver_param)
                #     update_model!(master, data)

                #     oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
                #     update_model!(oracle, data)

                #     root_seq_type = BendersSeq
                #     root_param = BendersSeqParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true
                #         )
                    
                #     # Create root node preprocessing with oracle
                #     root_preprocessing = NoRootNodePreprocessing()
                    
                #     # Create callbacks
                #     lazy_callback = LazyCallback(oracle=oracle)
                #     user_callback = NoUserCallback()
                    
                #     # Create BnB parameter
                #     callback_param = BendersBnBParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true
                #         )
                    
                #     # Create BendersBnB environment
                #     env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param=callback_param)
                    
                #     # Solve
                #     obj_value, elapsed_time = solve!(env)
                    
                #     @test env.termination_status == Optimal()
                #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                # end
                # @testset "Seq root preprocessing" begin        
                #     @info "solving p$i - CFLKnapsack oracle - seq..."
                #     master = Master(data; solver_param = master_solver_param)
                #     update_model!(master, data)

                #     oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
                #     update_model!(oracle, data)

                #     root_seq_type = BendersSeq
                #     root_param = BendersSeqParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true
                #         )
                    
                #     # Create root node preprocessing with oracle
                #     root_preprocessing = RootNodePreprocessing(oracle, root_seq_type, root_param)
                    
                #     # Create callbacks
                #     lazy_callback = LazyCallback(oracle=oracle)
                #     user_callback = UserCallback(params=UserCallbackParam(frequency=250), oracle=oracle)
                    
                #     # Create BnB parameter
                #     callback_param = BendersBnBParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true
                #         )
                    
                #     # Create BendersBnB environment
                #     env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param=callback_param)
                    
                #     # Solve
                #     obj_value, elapsed_time = solve!(env)
                    
                #     @test env.termination_status == Optimal()
                #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                # end
                # @testset "SeqInOut root preprocessing" begin        
                #     @info "solving p$i - CFLKnapsack oracle - seqinout..."
                #     master = Master(data; solver_param = master_solver_param)
                #     update_model!(master, data)

                #     oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
                #     update_model!(oracle, data)

                #     root_seq_type = BendersSeqInOut
                #     root_param = BendersSeqInOutParam(;
                #             time_limit = 100.0,
                #             gap_tolerance = 1e-6,
                #             stabilizing_x = ones(data.dim_x),
                #             α = 0.9,
                #             λ = 0.1,
                #             verbose = true
                #         )
                    
                #     # Create root node preprocessing with oracle
                #     root_preprocessing = RootNodePreprocessing(oracle, root_seq_type, root_param)
                    
                #     # Create callbacks
                #     lazy_callback = LazyCallback(oracle=oracle)
                #     user_callback = UserCallback(params=UserCallbackParam(frequency=250), oracle=oracle)
                    
                #     # Create BnB parameter
                #     callback_param = BendersBnBParam(;
                #             time_limit = 200.0,
                #             gap_tolerance = 1e-6,
                #             verbose = true
                #         )
                    
                #     # Create BendersBnB environment
                #     env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param=callback_param)
                    
                #     # Solve
                #     obj_value, elapsed_time = solve!(env)
                    
                #     @test env.termination_status == Optimal()
                #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                # end
            end

            # @testset "CFLKnapsack disjunctive oracle" begin
            #     @info "solving p$i - CFLKnapsack disjunctive oracle..."
            #     for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf], disjunctive_cut_append_rule in [NoDisjunctiveCuts(); AllDisjunctiveCuts(); DisjunctiveCutsSmallerIndices()]
            #         @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; p $p; dcut_append $disjunctive_cut_append_rule" begin
            #             master = Master(data; solver_param = master_solver_param)
            #             update_model!(master, data)

            #             lazy_oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
            #             update_model!(lazy_oracle, data)

            #             typical_oracles = [CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param); CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)] # for kappa & nu
            #             for k=1:2
            #                 update_model!(typical_oracles[k], data)
            #             end
            #             disjunctive_oracle = DisjunctiveOracle(data, typical_oracles; 
            #                             solver_param = dcglp_solver_param,
            #                             param = dcglp_param) 
            #             oracle_param = DisjunctiveOracleParam(norm = LpNorm(p), 
            #                 split_index_selection_rule = RandomFractional(),
            #                 disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
            #                 strengthened=strengthened, 
            #                 add_benders_cuts_to_master=add_benders_cuts_to_master, 
            #                 fraction_of_benders_cuts_to_master = 0.5, 
            #                 reuse_dcglp=reuse_dcglp)
            #             set_parameter!(disjunctive_oracle, oracle_param)
            #             update_model!(disjunctive_oracle, data)

            #             root_seq_type = BendersSeqInOut
            #             root_param = BendersSeqInOutParam(;
            #                 time_limit = 100.0,
            #                 gap_tolerance = 1e-6,
            #                 stabilizing_x = ones(data.dim_x),
            #                 α = 0.9,
            #                 λ = 0.1,
            #                 verbose = true
            #             )
                    
            #             # Create root node preprocessing with oracle
            #             root_preprocessing = RootNodePreprocessing(lazy_oracle, root_seq_type, root_param)
                        
            #             # Create callbacks
            #             lazy_callback = LazyCallback(oracle=lazy_oracle)
            #             user_callback = UserCallback(params=UserCallbackParam(frequency=1250), oracle=disjunctive_oracle)
            #             callback_param = BendersBnBParam(;
            #                 time_limit = 200.0,
            #                 gap_tolerance = 1e-6,
            #                 verbose = true
            #             )
            #             env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param=callback_param)
            #             log = solve!(env)
            #             @test env.termination_status == Optimal()
            #             @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            #         end
            #     end
            # end
        end
    end
end

