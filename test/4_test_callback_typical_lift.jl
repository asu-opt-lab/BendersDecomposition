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
import Random
global_logger(ConsoleLogger(stderr, Logging.Debug))
# Random.seed!(1218)

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
        halt_limit = 250, 
        iter_limit = 250,
        verbose = true
    )
    
    # Common solver parameters
    common_params = Dict(
        "solver" => "CPLEX", 
        "CPX_PARAM_EPINT" => 1e-9, 
        "CPX_PARAM_EPRHS" => 1e-9,
        "CPX_PARAM_EPGAP" => 1e-9
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

# Define function to setup and run standard test configurations
function run_standard_test(data, oracle, root_preproc_type, test_info)
    # Unpack parameters
    mip_opt_val, master_solver_param, typical_oracle_solver_param = test_info
    
    # Setup master problem
    master = Master(data; solver_param = master_solver_param)
    update_model!(master, data)
    
    # Setup root preprocessing based on type
    if root_preproc_type == :none
        root_preprocessing = NoRootNodePreprocessing()
        user_callback = NoUserCallback()
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
        root_preprocessing = RootNodePreprocessing(oracle, root_seq_type, root_param)
        user_callback = UserCallback(oracle; params=UserCallbackParam(frequency=250))
    end
    
    # Create callbacks
    lazy_callback = LazyCallback(oracle)
    
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
    
    # Test results
    @test env.termination_status == Optimal()
    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
    
    return env
end

# Setup and run disjunctive oracle test
function run_disjunctive_test(data, test_info, strengthened, add_benders_cuts_to_master, reuse_dcglp, p, disjunctive_cut_append_rule)
    # Unpack parameters
    mip_opt_val, master_solver_param, typical_oracle_solver_param, dcglp_solver_param, dcglp_param = test_info
    
    # Setup master problem
    master = Master(data; solver_param = master_solver_param)
    update_model!(master, data)

    # Create lazy oracle
    lazy_oracle = UFLKnapsackOracle(data)

    # Create typical oracles for kappa & nu
    typical_oracles = [
        UFLKnapsackOracle(data), 
        UFLKnapsackOracle(data)
    ]
    
    # Setup disjunctive oracle
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

    # Setup preprocessing
    root_seq_type = BendersSeq
    root_param = BendersSeqParam(;
        time_limit = 200.0,
        gap_tolerance = 1e-6,
        verbose = true
    )

    # Create root node preprocessing with oracle
    root_preprocessing = RootNodePreprocessing(lazy_oracle, root_seq_type, root_param)
    
    # Create callbacks
    lazy_callback = LazyCallback(lazy_oracle)
    user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=500))
    
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
    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
    
    return env
end

# Run test suite for a specific oracle type
function run_oracle_tests(data, oracle, test_info, root_preproc_types = [:none, :seq, :seqinout], description = "")
    for root_type in root_preproc_types
        root_type_str = string(root_type)
        @testset "$description - $root_type_str root preprocessing" begin
            @info "solving $(description) - $root_type_str root preprocessing..."
            # run_standard_test(data, oracle, root_type, test_info)
            run_disjunctive_test(data, test_info, true, false, false, Inf, NoDisjunctiveCuts())
            # data, test_info, strengthened, add_benders_cuts_to_master, reuse_dcglp, p, disjunctive_cut_append_rule
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
    instances = 4:4  # For quick testing
    
    for i in instances
        @testset "Instance: p$i" begin
            @info "Testing UFLP instance $i"
            
            # Load problem data
            # problem = read_uflp_benchmark_data("p49")
            problem = read_Simple_data("ga250a-2")
            
            # Get standard parameters
            benders_param, dcglp_param, mip_solver_param, master_solver_param, 
            typical_oracle_solver_param, dcglp_solver_param = get_standard_params()
            
            # Create data object for regular cuts
            # data = create_data(problem)
            
            # # Solve MIP for reference
            # mip_opt_val = solve_reference_mip(data, mip_solver_param)
            
            # Standard test info
            test_info = (0, master_solver_param, typical_oracle_solver_param, dcglp_solver_param, dcglp_param) # 0: arbitrary value for mip_opt_val
            
            # Test classical oracle
            # @testset "Classic oracle" begin
            #     oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param)
            #     update_model!(oracle, data)
            #     run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "classical oracle")
            # end
            
            # Create data object for knapsack cuts
            data = create_data(problem, problem.n_customers, ones(problem.n_customers))
            
            # Test fat knapsack oracle
            @testset "Fat knapsack oracle" begin
                oracle = UFLKnapsackOracle(data)
                set_parameter!(oracle, "add_only_violated_cuts", true)
                run_oracle_tests(data, oracle, test_info, [:seq], "fat knapsack oracle")
            end
            
            # Test slim knapsack oracle
            # @testset "Slim knapsack oracle" begin
            #     oracle = UFLKnapsackOracle(data)
            #     set_parameter!(oracle, "add_only_violated_cuts", false)
            #     set_parameter!(oracle, "slim", true)
            #     run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "slim knapsack oracle")
            # end
        end
    end
end