using Test
using JuMP
using CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition
import BendersDecomposition: generate_cuts


function get_standard_params()
    # Algorithm parameters
    benders_param = BendersBnBParam(;
        time_limit = 200.0,
        gap_tolerance = 1e-6,
        verbose = false
    )
    
    dcglp_param = DcglpParam(;
        time_limit = 1000.0, 
        gap_tolerance = 1e-3, 
        halt_limit = 3, 
        iter_limit = 250,
        verbose = false
    )
    
    # Common solver parameters
    common_params = Dict(
        "solver" => "CPLEX", 
        "CPX_PARAM_EPINT" => 1e-9, 
        "CPX_PARAM_EPRHS" => 1e-9,
        "CPXPARAM_Threads" => 4
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
                verbose = false
            )
        elseif root_preproc_type == :seqinout
            root_seq_type = BendersSeqInOut
            root_param = BendersSeqInOutParam(;
                time_limit = root_preproc_type == :seqinout ? 100.0 : 200.0,
                gap_tolerance = 1e-6,
                stabilizing_x = ones(data.dim_x),
                α = 0.9,
                λ = 0.1,
                verbose = false
            )
        end
        
        # Create root node preprocessing with oracle
        root_preprocessing = RootNodePreprocessing(oracle, root_seq_type, root_param)
        user_callback = UserCallback(oracle; params=UserCallbackParam(frequency=10))
    end
    
    # Create callbacks
    lazy_callback = LazyCallback(oracle)
    
    # Create BnB parameter
    callback_param = BendersBnBParam(;
        time_limit = 200.0,
        gap_tolerance = 1e-6,
        verbose = false
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

# Run test suite for a specific oracle type
function run_oracle_tests(data, oracle, test_info, root_preproc_types = [:none, :seq, :seqinout], description = "")
    for root_type in root_preproc_types
        root_type_str = string(root_type)
        @testset "$description - $root_type_str root preprocessing" begin
            @info "solving $(description) - $root_type_str root preprocessing..."
            run_standard_test(data, oracle, root_type, test_info)
        end
    end
end

@testset "Callback Typical Tests" begin
    @info "Running Callback Typical Tests"
    include("ufl.jl")
    include("cfl.jl")
    include("scfl.jl")
    @info "Callback Typical Tests completed"
end