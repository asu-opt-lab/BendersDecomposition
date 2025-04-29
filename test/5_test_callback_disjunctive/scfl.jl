include("$(dirname(dirname(@__DIR__)))/example/scflp/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/scflp/model.jl")

function create_scflp_oracle(data, oracle_type, oracle_solver_param)
    oracle = SeparableOracle(data, oracle_type(), data.problem.n_scenarios; solver_param = oracle_solver_param)
    for j=1:oracle.N
        update_model!(oracle.oracles[j], data, j)
    end
    return oracle
end

@testset verbose = true "SCFLP Callback Benders Tests" begin
    # Specify instances to test
    instances = 1:1  # For quick testing
    
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
            
            @testset "Classic oracle" begin
                # Create separable oracle
                typical_oracle_kappa = create_scflp_oracle(data, ClassicalOracle, typical_oracle_solver_param)
                typical_oracle_nu = create_scflp_oracle(data, ClassicalOracle, typical_oracle_solver_param)
                typical_oracles = [typical_oracle_kappa; typical_oracle_nu]
                
                # Create lazy oracle
                lazy_oracle = create_scflp_oracle(data, ClassicalOracle, typical_oracle_solver_param)

                for strengthened in [true, false], 
                    add_benders_cuts_to_master in [true, false], 
                    reuse_dcglp in [true, false], 
                    p in [1.0, Inf], 
                    disjunctive_cut_append_rule in [NoDisjunctiveCuts(), AllDisjunctiveCuts(), DisjunctiveCutsSmallerIndices()]
                    
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; p $p; dcut_append $disjunctive_cut_append_rule" begin
                        @info "solving SCFLP p$i - disjunctive oracle/classical - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule"

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
                            fraction_of_benders_cuts_to_master = 1.0, 
                            reuse_dcglp = reuse_dcglp
                        )
                        set_parameter!(disjunctive_oracle, oracle_param)
                        update_model!(disjunctive_oracle, data)

                        run_disjunctive_oracle_tests(data, mip_opt_val, lazy_oracle, disjunctive_oracle, [:none, :seq, :seqinout], "SCFLP Classical oracle")
                    end
                end
            end

            @testset "Knapsack oracle" begin
                # Create separable oracle
                typical_oracle_kappa = create_scflp_oracle(data, CFLKnapsackOracle, typical_oracle_solver_param)
                typical_oracle_nu = create_scflp_oracle(data, CFLKnapsackOracle, typical_oracle_solver_param)
                typical_oracles = [typical_oracle_kappa; typical_oracle_nu]
                
                # Create lazy oracle
                lazy_oracle = create_scflp_oracle(data, CFLKnapsackOracle, typical_oracle_solver_param)

                for strengthened in [true, false], 
                    add_benders_cuts_to_master in [true, false], 
                    reuse_dcglp in [true, false], 
                    p in [1.0, Inf], 
                    disjunctive_cut_append_rule in [NoDisjunctiveCuts(), AllDisjunctiveCuts(), DisjunctiveCutsSmallerIndices()]
                    
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; p $p; dcut_append $disjunctive_cut_append_rule" begin
                        @info "solving SCFLP p$i - disjunctive oracle/knapsack - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule"

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
                            fraction_of_benders_cuts_to_master = 1.0, 
                            reuse_dcglp = reuse_dcglp
                        )
                        set_parameter!(disjunctive_oracle, oracle_param)
                        update_model!(disjunctive_oracle, data)

                        run_disjunctive_oracle_tests(data, mip_opt_val, lazy_oracle, disjunctive_oracle, [:none, :seq, :seqinout], "SCFLP Knapsack oracle")
                    end
                end
            end
        end
    end
end


            
            
            
            
            




