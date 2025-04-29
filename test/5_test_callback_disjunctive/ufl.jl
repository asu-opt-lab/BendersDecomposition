include("$(dirname(dirname(@__DIR__)))/example/uflp/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/uflp/oracle.jl")
include("$(dirname(dirname(@__DIR__)))/example/uflp/model.jl")

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
            
            # keep this for future reference
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