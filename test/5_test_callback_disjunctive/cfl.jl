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
                @info "solving CFLP p$i - classical oracle..."
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
                    add_benders_cuts_to_master in [true, false, 2], 
                    reuse_dcglp in [true, false], 
                    lift in [true, false],
                    p in [1.0, Inf], 
                    disjunctive_cut_append_rule in [NoDisjunctiveCuts(), AllDisjunctiveCuts(), DisjunctiveCutsSmallerIndices()]
                    
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; lift $lift; p $p; dcut_append $disjunctive_cut_append_rule" begin
                        @info "solving CFLP p$i - disjunctive oracle/classical - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule"
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
                            reuse_dcglp = reuse_dcglp,
                            lift = lift
                        )
                        set_parameter!(disjunctive_oracle, oracle_param)
                        update_model!(disjunctive_oracle, data)
                        
                        run_disjunctive_oracle_tests(data, mip_opt_val, lazy_oracle, disjunctive_oracle, [:none, :seq, :seqinout], "CFLP Classical oracle")
                    end
                end
            end
            
            # Test CFLKnapsack oracle
            @testset "CFLKnapsack oracle" begin
                lazy_oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
                update_model!(lazy_oracle, data)
                typical_oracles = [
                    CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param), 
                    CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
                ]
                for k=1:2
                    update_model!(typical_oracles[k], data)
                end
                
                # Test various parameter combinations
                for strengthened in [true, false], 
                    add_benders_cuts_to_master in [true, false, 2],
                    reuse_dcglp in [true, false],
                    lift in [true, false],
                    p in [1.0, Inf],
                    disjunctive_cut_append_rule in [NoDisjunctiveCuts(), AllDisjunctiveCuts(), DisjunctiveCutsSmallerIndices()]
                    
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; lift $lift; p $p; dcut_append $disjunctive_cut_append_rule" begin
                        @info "solving CFLP p$i - disjunctive oracle/knapsack - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule"
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
                            reuse_dcglp = reuse_dcglp,
                            lift = lift
                        )
                        set_parameter!(disjunctive_oracle, oracle_param)
                        update_model!(disjunctive_oracle, data)

                        run_disjunctive_oracle_tests(data, mip_opt_val, lazy_oracle, disjunctive_oracle, [:none, :seq, :seqinout], "CFLKnapsack oracle")
                    end
                end
            end
        end
    end

    instances = 1:5
    for i in instances
        @testset "Instance: p$i" begin
            @info "Testing CFLP hard instance $i"
            
            # Load problem data
            problem = read_GK_data("f100-c100-r5-$i")
            
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
                @info "solving CFLP p$i - classical oracle..."
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
                    add_benders_cuts_to_master in [true, false, 2], 
                    reuse_dcglp in [true, false], 
                    lift in [true, false],
                    p in [1.0, Inf], 
                    disjunctive_cut_append_rule in [NoDisjunctiveCuts(), AllDisjunctiveCuts(), DisjunctiveCutsSmallerIndices()]
                    
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; lift $lift p $p; dcut_append $disjunctive_cut_append_rule" begin
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
                        
                        run_disjunctive_oracle_tests(data, mip_opt_val, lazy_oracle, disjunctive_oracle, [:none, :seq, :seqinout], "CFL Classical oracle")
                    end
                end
            end
            
            # Test CFLKnapsack oracle
            @testset "CFLKnapsack oracle" begin
                lazy_oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
                update_model!(lazy_oracle, data)
                typical_oracles = [
                    CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param), 
                    CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
                ]
                for k=1:2
                    update_model!(typical_oracles[k], data)
                end
                
                # Test various parameter combinations
                for strengthened in [true, false], 
                    add_benders_cuts_to_master in [true, false, 2],
                    reuse_dcglp in [true, false],
                    lift in [true, false],
                    p in [1.0, Inf],
                    disjunctive_cut_append_rule in [NoDisjunctiveCuts(), AllDisjunctiveCuts(), DisjunctiveCutsSmallerIndices()]
                    
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; lift $lift p $p; dcut_append $disjunctive_cut_append_rule" begin
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
                            reuse_dcglp = reuse_dcglp,
                            lift = lift
                        )
                        set_parameter!(disjunctive_oracle, oracle_param)
                        update_model!(disjunctive_oracle, data)

                        run_disjunctive_oracle_tests(data, mip_opt_val, lazy_oracle, disjunctive_oracle, [:none, :seq, :seqinout], "CFLKnapsack oracle")
                    end
                end
            end
        end
    end
end