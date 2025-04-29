# Include file dependencies
include("$(dirname(dirname(@__DIR__)))/example/uflp/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/uflp/oracle.jl")
include("$(dirname(dirname(@__DIR__)))/example/uflp/model.jl")

@testset verbose = true "UFLP Callback Benders Tests" begin
    # Specify instances to test
    instances = setdiff(1:71, [67])  # For quick testing
    
    for i in instances
        @testset "Instance: p$i" begin
            @info "Testing UFLP instance $i"
            
            # Load problem data
            problem = read_uflp_benchmark_data("p49")
            
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
                oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param)
                update_model!(oracle, data)
                run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "classical oracle")
            end
            
            # Create data object for knapsack cuts
            data = create_data(problem, problem.n_customers, ones(problem.n_customers))
            
            # Test fat knapsack oracle
            @testset "Fat knapsack oracle" begin
                oracle = UFLKnapsackOracle(data)
                set_parameter!(oracle, "add_only_violated_cuts", true)
                run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "fat knapsack oracle")
            end
            
            # keep this for future reference
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