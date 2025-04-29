include("$(dirname(dirname(@__DIR__)))/example/cflp/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/cflp/oracle.jl")
include("$(dirname(dirname(@__DIR__)))/example/cflp/model.jl")

@testset verbose = true "CFLP Callback Benders Tests" begin
    # Specify instances to test
    instances = setdiff(1:71, [67])  # For quick testing
    
    for i in instances
        @testset "Instance: p$i" begin
            @info "Testing CFLP easy instance $i"
            
            # Load problem data
            problem = read_cflp_benchmark_data("p$i")
            
            # Get standard parameters
            benders_param, dcglp_param, mip_solver_param, master_solver_param, 
            typical_oracle_solver_param, dcglp_solver_param = get_standard_params()
            
            # Create data object
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
            
            # Test CFLKnapsack oracle
            @testset "CFLKnapsack oracle" begin
                oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
                update_model!(oracle, data)
                run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "CFLKnapsack oracle")
            end
        end
    end

    # instances = 1:5
    # for i in instances
    #     @testset "Instance: p$i" begin
    #         @info "Testing CFLP hard instance $i"
            
    #         # Load problem data
    #         problem = read_GK_data("f100-c100-r5-$i")
            
    #         # Get standard parameters
    #         benders_param, dcglp_param, mip_solver_param, master_solver_param, 
    #         typical_oracle_solver_param, dcglp_solver_param = get_standard_params()
            
    #         # Create data object
    #         data = create_data(problem)
            
    #         # Solve MIP for reference
    #         mip_opt_val = solve_reference_mip(data, mip_solver_param)
            
    #         # Standard test info
    #         test_info = (mip_opt_val, master_solver_param, typical_oracle_solver_param)
            
    #         # Test classical oracle
    #         @testset "Classic oracle" begin
    #             oracle = ClassicalOracle(data; solver_param = typical_oracle_solver_param)
    #             update_model!(oracle, data)
    #             run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "classical oracle")
    #         end
            
    #         # Test CFLKnapsack oracle
    #         @testset "CFLKnapsack oracle" begin
    #             oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
    #             update_model!(oracle, data)
    #             run_oracle_tests(data, oracle, test_info, [:none, :seq, :seqinout], "CFLKnapsack oracle")
    #         end
    #     end
    # end
end