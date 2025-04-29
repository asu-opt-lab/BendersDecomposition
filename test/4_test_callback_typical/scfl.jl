# Include SCFLP model files
include("$(dirname(dirname(@__DIR__)))/example/scflp/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/scflp/model.jl")

# Create and initialize separable oracle for SCFLP
function create_scflp_oracle(data, oracle_type, oracle_solver_param)
    oracle = SeparableOracle(data, oracle_type(), data.problem.n_scenarios; solver_param = oracle_solver_param)
    for j=1:oracle.N
        update_model!(oracle.oracles[j], data, j)
    end
    return oracle
end

@testset verbose = true "SCFLP Sequential Benders Tests" begin
    # Specify instances to test
    instances = 1:5  # For quick testing
    
    for i in instances
        @testset "Instance: f25-c50-s64-r10-$i" begin
            @info "Testing SCFLP instance $i"
            
            # Load problem data
            problem = read_stochastic_capacited_facility_location_problem("f25-c50-s64-r10-$i")
            
            # Initialize data object
            dim_x = problem.n_facilities
            dim_t = problem.n_scenarios
            c_x = problem.fixed_costs
            c_t = fill(1/problem.n_scenarios, problem.n_scenarios)
            data = Data(dim_x, dim_t, problem, c_x, c_t)
            
            # Get standard parameters
            benders_param, dcglp_param, mip_solver_param, master_solver_param, 
            typical_oracle_solver_param, dcglp_solver_param = get_standard_params()
            
            # Solve MIP for reference
            mip_opt_val = solve_reference_mip(data, mip_solver_param)
            
            # Standard test info
            test_info = (mip_opt_val, master_solver_param, typical_oracle_solver_param)
            
            # Test classical oracle
            @testset "Classic oracle" begin
                oracle = create_scflp_oracle(data, ClassicalOracle, typical_oracle_solver_param)
                run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "classical oracle")
            end 
            
            # Test knapsack oracle
            @testset "Knapsack oracle" begin
                oracle = create_scflp_oracle(data, CFLKnapsackOracle, typical_oracle_solver_param)
                run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "knapsack oracle")
            end
        end
    end
end