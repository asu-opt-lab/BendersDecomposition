# To-Do: 
# 1. Need to be able to change the setting for SeqInOut (e.g., stabilizing point)
# 2. Assign attributes to the structure, not to JuMP Model (e.g., settings for CFLKnapsackOracle like slim, add_only_violated_cuts)
# Done:
# - Slim cut should be averaged, instead of summation, for numerical stability

using Test
using JuMP
using CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition
import BendersDecomposition: generate_cuts

# -----------------------------------------------------------------------------
# Common test utilities and parameter settings
# -----------------------------------------------------------------------------

# Get standard solver and algorithm parameters
function get_standard_params()
    # Algorithm parameters
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
    
    # Common solver parameters
    common_params = Dict(
        "solver" => "CPLEX", 
        "CPX_PARAM_EPINT" => 1e-9, 
        "CPX_PARAM_EPRHS" => 1e-9
    )
    
    # Oracle-specific parameters
    oracle_params = Dict(
        "solver" => "CPLEX", 
        "CPX_PARAM_EPRHS" => 1e-9, 
        "CPX_PARAM_NUMERICALEMPHASIS" => 1, 
        "CPX_PARAM_EPOPT" => 1e-9
    )
    
    return benders_param, dcglp_param, common_params, common_params, oracle_params, oracle_params
end

# Create data object for a given problem
function create_data(problem, dim_t = 1, c_t = [1])
    dim_x = problem.n_facilities
    c_x = problem.fixed_costs
    
    data = Data(dim_x, dim_t, problem, c_x, c_t)
    @assert dim_x == length(data.c_x)
    @assert dim_t == length(data.c_t)
    
    return data
end

# Solve MIP for reference
function solve_reference_mip(data, mip_solver_param)
    mip = Mip(data)
    assign_attributes!(mip.model, mip_solver_param)
    update_model!(mip, data)
    optimize!(mip.model)
    @assert termination_status(mip.model) == OPTIMAL
    return objective_value(mip.model)
end

# Setup and run disjunctive oracle test
function run_disjunctive_test(data, lazy_oracle, disjunctive_oracle, root_preproc_type)
    
    # Setup master problem
    master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9)
    master = Master(data; solver_param = master_solver_param)
    update_model!(master, data)

    # Setup preprocessing
    if root_preproc_type == :none
        root_preprocessing = NoRootNodePreprocessing()
        user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=250))
    else
        if root_preproc_type == :seq
            root_seq_type = BendersSeq
            root_param = BendersSeqParam(;
                time_limit = 200.0,
                gap_tolerance = 1e-6,
                verbose = true
            )
        elseif root_preproc_type == :seqinout
            root_seq_type = BendersSeqInOut
            root_param = BendersSeqInOutParam(;
                time_limit = root_preproc_type == :seqinout ? 100.0 : 200.0,
                gap_tolerance = 1e-6,
                stabilizing_x = ones(data.dim_x),
                α = 0.9,
                λ = 0.1,
                verbose = true
            )
        end
        
        # Create root node preprocessing with oracle
        root_preprocessing = RootNodePreprocessing(lazy_oracle, root_seq_type, root_param)
        user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=250))
    end
    
    lazy_callback = LazyCallback(lazy_oracle)
    # Create BnB parameter
    callback_param = BendersBnBParam(;
        time_limit = 200.0,
        gap_tolerance = 1e-6,
        verbose = true
    )
    
    # Create BendersBnB environment
    env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param=callback_param)
    
    # Solve
    log = solve!(env)
    
    # Test results
    @test env.termination_status == Optimal()
    
    return env
end

# Run test suite for a specific oracle type
function run_disjunctive_oracle_tests(data, mip_opt_val, lazy_oracle, disjunctive_oracle, root_preproc_types = [:none, :seq, :seqinout], description = "")
    for root_type in root_preproc_types
        root_type_str = string(root_type)
        @testset "$description - $root_type_str root preprocessing" begin
            @info "solving $(description) - $root_type_str root preprocessing..."
            env = run_disjunctive_test(data, lazy_oracle, disjunctive_oracle, root_type)
            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
        end
    end
end

# -----------------------------------------------------------------------------
# UFLP Tests
# -----------------------------------------------------------------------------

# Include file dependencies
include("$(dirname(@__DIR__))/example/uflp/data_reader.jl")
include("$(dirname(@__DIR__))/example/uflp/oracle.jl")
include("$(dirname(@__DIR__))/example/uflp/model.jl")

@testset verbose = true "UFLP Callback Benders Tests" begin
    # Specify instances to test
    instances = setdiff(1:71, [67])  # For quick testing
    # instances = 25:25

    for i in instances
        @testset "Instance: p$i" begin
            @info "Testing UFLP instance $i"
            
            # Load problem data
            problem = read_uflp_benchmark_data("p$i")
            
            # Get standard parameters
            benders_param, dcglp_param, mip_solver_param, master_solver_param, 
            typical_oracle_solver_param, dcglp_solver_param = get_standard_params()
            
            # Create data object for regular cuts
            data = create_data(problem)
            
            # Solve MIP for reference
            mip_opt_val = solve_reference_mip(data, mip_solver_param)
            
            # Standard test info
            test_info = (mip_opt_val, master_solver_param, typical_oracle_solver_param)
            
            # Test classical oracle
            @testset "Classic oracle" begin
                lazy_oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param)
                update_model!(lazy_oracle, data)
                typical_oracles = [
                    ClassicalOracle(data; solver_param = typical_oracle_solver_param), 
                    ClassicalOracle(data; solver_param = typical_oracle_solver_param)
                ]
                for k=1:2
                    update_model!(typical_oracles[k], data)
                end
                
                # Test various parameter combinations
                for strengthened in [true, false], 
                    add_benders_cuts_to_master in [true, false], 
                    reuse_dcglp in [true, false], 
                    p in [1.0, Inf], 
                    disjunctive_cut_append_rule in [NoDisjunctiveCuts(), AllDisjunctiveCuts(), DisjunctiveCutsSmallerIndices()]
                    
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; p $p; dcut_append $disjunctive_cut_append_rule" begin
                        disjunctive_oracle = DisjunctiveOracle(data, typical_oracles; 
                            solver_param = dcglp_solver_param,
                            param = dcglp_param
                        ) 

                        # Set oracle parameters
                        oracle_param = DisjunctiveOracleParam(
                            norm = LpNorm(p), 
                            split_index_selection_rule = RandomFractional(),
                            disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
                            strengthened = strengthened, 
                            add_benders_cuts_to_master = add_benders_cuts_to_master, 
                            fraction_of_benders_cuts_to_master = 0.5, 
                            reuse_dcglp = reuse_dcglp
                        )
                        set_parameter!(disjunctive_oracle, oracle_param)
                        update_model!(disjunctive_oracle, data)
                        
                        run_disjunctive_oracle_tests(data, mip_opt_val, lazy_oracle, disjunctive_oracle, [:none, :seq, :seqinout], "UFLP")
                    end
                end
            end
            
            # Create data object for knapsack cuts
            data = create_data(problem, problem.n_customers, ones(problem.n_customers))
            
            # Test fat knapsack oracle
            @testset "Fat knapsack oracle" begin
                lazy_oracle = UFLKnapsackOracle(data)
                set_parameter!(lazy_oracle, "add_only_violated_cuts", true)
                typical_oracles = [
                    UFLKnapsackOracle(data), 
                    UFLKnapsackOracle(data)
                ]
                for k=1:2
                    set_parameter!(typical_oracles[k], "add_only_violated_cuts", true)
                end

                for strengthened in [true, false], 
                    add_benders_cuts_to_master in [true, false], 
                    reuse_dcglp in [true, false], 
                    p in [1.0, Inf], 
                    disjunctive_cut_append_rule in [NoDisjunctiveCuts(), AllDisjunctiveCuts(), DisjunctiveCutsSmallerIndices()]
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; p $p; dcut_append $disjunctive_cut_append_rule" begin
                        disjunctive_oracle = DisjunctiveOracle(data, typical_oracles; 
                            solver_param = dcglp_solver_param,
                            param = dcglp_param
                        ) 

                        # Set oracle parameters
                        oracle_param = DisjunctiveOracleParam(
                            norm = LpNorm(p), 
                            split_index_selection_rule = RandomFractional(),
                            disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
                            strengthened = strengthened, 
                            add_benders_cuts_to_master = add_benders_cuts_to_master, 
                            fraction_of_benders_cuts_to_master = 0.5, 
                            reuse_dcglp = reuse_dcglp
                        )
                        set_parameter!(disjunctive_oracle, oracle_param)
                        update_model!(disjunctive_oracle, data)
                        
                        run_disjunctive_oracle_tests(data, mip_opt_val, lazy_oracle, disjunctive_oracle, [:none, :seq, :seqinout], "UFLP")
                    end
                end
            end
            
            # Test slim knapsack oracle
            # @testset "Slim knapsack oracle" begin
            #     lazy_oracle = UFLKnapsackOracle(data)
            #     set_parameter!(lazy_oracle, "add_only_violated_cuts", false)
            #     set_parameter!(lazy_oracle, "slim", true)
                
            #     typical_oracles = [
            #         UFLKnapsackOracle(data), 
            #         UFLKnapsackOracle(data)
            #     ]
            #     for k=1:2
            #         set_parameter!(typical_oracles[k], "add_only_violated_cuts", false)
            #         set_parameter!(typical_oracles[k], "slim", true)
            #     end

            #     for strengthened in [true, false], 
            #         add_benders_cuts_to_master in [true, false], 
            #         reuse_dcglp in [true, false], 
            #         p in [1.0, Inf], 
            #         disjunctive_cut_append_rule in [NoDisjunctiveCuts(), AllDisjunctiveCuts(), DisjunctiveCutsSmallerIndices()]
            #         @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; p $p; dcut_append $disjunctive_cut_append_rule" begin
            #             disjunctive_oracle = DisjunctiveOracle(data, typical_oracles; 
            #                 solver_param = dcglp_solver_param,
            #                 param = dcglp_param
            #             ) 

            #             oracle_param = DisjunctiveOracleParam(
            #                 norm = LpNorm(p), 
            #                 split_index_selection_rule = RandomFractional(),
            #                 disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
            #                 strengthened = strengthened, 
            #                 add_benders_cuts_to_master = add_benders_cuts_to_master, 
            #                 fraction_of_benders_cuts_to_master = 0.5, 
            #                 reuse_dcglp = reuse_dcglp
            #             )
            #             set_parameter!(disjunctive_oracle, oracle_param)
            #             update_model!(disjunctive_oracle, data)

            #             run_disjunctive_oracle_tests(data, mip_opt_val, lazy_oracle, disjunctive_oracle, [:none, :seq, :seqinout], "UFLP")
            #         end
            #     end
            # end
        end
    end
end

# -----------------------------------------------------------------------------
# CFLP Tests
# -----------------------------------------------------------------------------

# # Include CFLP model files 
# include("$(dirname(@__DIR__))/example/cflp/data_reader.jl")
# include("$(dirname(@__DIR__))/example/cflp/oracle.jl")
# include("$(dirname(@__DIR__))/example/cflp/model.jl")

# @testset verbose = true "CFLP Callback Benders Tests" begin
#     # Specify instances to test
#     instances = setdiff(1:71, [67])  # For quick testing
    
#     for i in instances
#         @testset "Instance: p$i" begin
#             @info "Testing CFLP easy instance $i"
            
#             # Load problem data
#             problem = read_cflp_benchmark_data("p$i")
            
#             # Get standard parameters
#             benders_param, dcglp_param, mip_solver_param, master_solver_param, 
#             typical_oracle_solver_param, dcglp_solver_param = get_standard_params()
            
#             # Create data object
#             data = create_data(problem)
            
#             # Solve MIP for reference
#             mip_opt_val = solve_reference_mip(data, mip_solver_param)
            
#             # Standard test info
#             test_info = (mip_opt_val, master_solver_param, typical_oracle_solver_param)
            
#             # Test classical oracle
#             @testset "Classic oracle" begin
#                 oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param)
#                 update_model!(oracle, data)
#                 run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "classical oracle")
#             end
            
#             # Test CFLKnapsack oracle
#             @testset "CFLKnapsack oracle" begin
#                 oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
#                 update_model!(oracle, data)
#                 run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "CFLKnapsack oracle")
#             end
            
#             Disjunctive oracle tests - commented out for brevity
#             @testset "CFLKnapsack disjunctive oracle" begin
#                 @info "solving p$i - CFLKnapsack disjunctive oracle..."
#                 disjunctive_test_info = (mip_opt_val, master_solver_param, typical_oracle_solver_param, dcglp_solver_param, dcglp_param)
                
#                 # Test various parameter combinations
#                 for strengthened in [true, false], 
#                     add_benders_cuts_to_master in [true, false], 
#                     reuse_dcglp in [true, false], 
#                     p in [1.0, Inf], 
#                     disjunctive_cut_append_rule in [NoDisjunctiveCuts(), AllDisjunctiveCuts(), DisjunctiveCutsSmallerIndices()]
                    
#                     @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; p $p; dcut_append $disjunctive_cut_append_rule" begin
#                         run_disjunctive_test(data, disjunctive_test_info, strengthened, add_benders_cuts_to_master, reuse_dcglp, p, disjunctive_cut_append_rule)
#                     end
#                 end
#             end
#         end
#     end

#     # instances = 1:5
#     # for i in instances
#     #     @testset "Instance: p$i" begin
#     #         @info "Testing CFLP hard instance $i"
            
#     #         # Load problem data
#     #         problem = read_GK_data("f100-c100-r5-$i")
            
#     #         # Get standard parameters
#     #         benders_param, dcglp_param, mip_solver_param, master_solver_param, 
#     #         typical_oracle_solver_param, dcglp_solver_param = get_standard_params()
            
#     #         # Create data object
#     #         data = create_data(problem)
            
#     #         # Solve MIP for reference
#     #         mip_opt_val = solve_reference_mip(data, mip_solver_param)
            
#     #         # Standard test info
#     #         test_info = (mip_opt_val, master_solver_param, typical_oracle_solver_param)
            
#     #         # Test classical oracle
#     #         @testset "Classic oracle" begin
#     #             oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param)
#     #             update_model!(oracle, data)
#     #             # run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "classical oracle")
#     #             run_oracle_tests(data, oracle, test_info, [:seqinout], "classical oracle")
#     #         end
            
#     #         # Test CFLKnapsack oracle
#     #         @testset "CFLKnapsack oracle" begin
#     #             oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
#     #             update_model!(oracle, data)
#     #             # run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "CFLKnapsack oracle")
#     #             run_oracle_tests(data, oracle, test_info, [:seqinout], "CFLKnapsack oracle")
#     #         end
            
#     #         # Disjunctive oracle tests - commented out for brevity
#     #         @testset "CFLKnapsack disjunctive oracle" begin
#     #             @info "solving p$i - CFLKnapsack disjunctive oracle..."
#     #             disjunctive_test_info = (mip_opt_val, master_solver_param, typical_oracle_solver_param, dcglp_solver_param, dcglp_param)
                
#     #             # Test various parameter combinations
#     #             for strengthened in [true, false], 
#     #                 add_benders_cuts_to_master in [true, false], 
#     #                 reuse_dcglp in [true, false], 
#     #                 p in [1.0, Inf], 
#     #                 disjunctive_cut_append_rule in [NoDisjunctiveCuts(), AllDisjunctiveCuts(), DisjunctiveCutsSmallerIndices()]
                    
#     #                 @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; p $p; dcut_append $disjunctive_cut_append_rule" begin
#     #                     run_disjunctive_test(data, disjunctive_test_info, strengthened, add_benders_cuts_to_master, reuse_dcglp, p, disjunctive_cut_append_rule)
#     #                 end
#     #             end
#     #         end
#     #     end
#     # end
# end



# -----------------------------------------------------------------------------
# SCFLP Tests
# -----------------------------------------------------------------------------

# Include SCFLP model files
# include("$(dirname(@__DIR__))/example/scflp/data_reader.jl")
# include("$(dirname(@__DIR__))/example/scflp/model.jl")

# # Create and initialize separable oracle for SCFLP
# function create_scflp_oracle(data, oracle_type, oracle_solver_param)
#     oracle = SeparableOracle(data, oracle_type(), data.problem.n_scenarios; solver_param = oracle_solver_param)
#     for j=1:oracle.N
#         update_model!(oracle.oracles[j], data, j)
#     end
#     return oracle
# end

# @testset verbose = true "SCFLP Sequential Benders Tests" begin
#     # Specify instances to test
#     instances = 1:5  # For quick testing
    
#     for i in instances
#         @testset "Instance: f25-c50-s64-r10-$i" begin
#             @info "Testing SCFLP instance $i"
            
#             # Load problem data
#             problem = read_stochastic_capacited_facility_location_problem("f25-c50-s64-r10-$i")
            
#             # Initialize data object
#             dim_x = problem.n_facilities
#             dim_t = problem.n_scenarios
#             c_x = problem.fixed_costs
#             c_t = fill(1/problem.n_scenarios, problem.n_scenarios)
#             data = Data(dim_x, dim_t, problem, c_x, c_t)
            
#             # Get standard parameters
#             benders_param, dcglp_param, mip_solver_param, master_solver_param, 
#             typical_oracle_solver_param, dcglp_solver_param = get_standard_params()
            
#             # Solve MIP for reference
#             mip_opt_val = solve_reference_mip(data, mip_solver_param)
            
#             # Standard test info
#             test_info = (mip_opt_val, master_solver_param, typical_oracle_solver_param)
            
#             # Test classical oracle
#             @testset "Classic oracle" begin
#                 oracle = create_scflp_oracle(data, ClassicalOracle, typical_oracle_solver_param)
#                 run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "classical oracle")
#             end 
            
#             # Test knapsack oracle
#             @testset "Knapsack oracle" begin
#                 oracle = create_scflp_oracle(data, CFLKnapsackOracle, typical_oracle_solver_param)
#                 run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "knapsack oracle")
#             end
#         end
#     end
# end