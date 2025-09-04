include("$(dirname(dirname(@__DIR__)))/example/uflp/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/uflp/oracle.jl")
include("$(dirname(dirname(@__DIR__)))/example/uflp/model.jl")

@testset verbose = true "UFLP Callback Disjunctive Benders Tests" begin
    # Specify instances to test
    instances = setdiff(1:71, [67])  # For quick testing
    # instances = 1:1

    for i in instances
        @testset "Instance: p$i" begin
            @info "Testing UFLP instance $i"
            
            # Load problem data
            problem = read_uflp_benchmark_data("p$i")
            
            # Get standard parameters
            # Get standard parameters
            benders_param = BendersBnBParam(;
                time_limit = 200.0,
                gap_tolerance = 1e-6,
                verbose = false
            )

            dcglp_param = DcglpParam(
                time_limit = 200.0,
                gap_tolerance = 1e-3,
                halt_limit = 3,
                iter_limit = 15,
                verbose = false
            )
            
            # Common solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPXPARAM_Threads" => 4)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 4)
            typical_oracle_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)
            dcglp_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, 
            "CPX_PARAM_EPOPT" => 1e-9) 
            
            # Create data object for regular cuts
            dim_x = problem.n_facilities
            dim_t = 1
            c_x = problem.fixed_costs
            c_t = [1]
            data = Data(dim_x, dim_t, problem, c_x, c_t)
            @assert dim_x == length(data.c_x)
            @assert dim_t == length(data.c_t)
            
            # Solve MIP for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)
            
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
                    lift = [true, false],
                    p in [1.0, Inf], 
                    disjunctive_cut_append_rule in [NoDisjunctiveCuts(), AllDisjunctiveCuts(), DisjunctiveCutsSmallerIndices()],
                    adjust_t_to_fx in [true; false]
                    
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; lift $lift; p $p; dcut_append $disjunctive_cut_append_rule; adjust_t_to_fx $adjust_t_to_fx" begin
                        @info "solving UFLP p$i - disjunctive oracle/classical - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"   
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
                            lift = lift,
                            adjust_t_to_fx = adjust_t_to_fx
                        )
                        set_parameter!(disjunctive_oracle, oracle_param)
                        update_model!(disjunctive_oracle, data)
                        
                        @testset "NoSeq" begin
                            @info "solving UFLP p$i - disjunctive oracle/classical/no seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"
                            master = Master(data; solver_param = master_solver_param)
                            update_model!(master, data)
                            root_preprocessing = NoRootNodePreprocessing()
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=10))
                            env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end
                        @testset "Seq" begin
                            @info "solving UFLP p$i - disjunctive oracle/classical/seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"
                            master = Master(data; solver_param = master_solver_param)
                            update_model!(master, data)
                            root_preprocessing = RootNodePreprocessing(lazy_oracle, BendersSeq, BendersSeqParam(;time_limit=200.0, gap_tolerance=1e-6, verbose=false))
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=10))
                            env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end
                        @testset "SeqInOut" begin
                            @info "solving UFLP p$i - disjunctive oracle/classical/seqinout - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"
                            master = Master(data; solver_param = master_solver_param)
                            update_model!(master, data)
                            root_preprocessing = RootNodePreprocessing(lazy_oracle, BendersSeqInOut, BendersSeqInOutParam(time_limit = 300.0, gap_tolerance = 1e-6, stabilizing_x = ones(data.dim_x), α = 0.9, λ = 0.1, verbose = false))
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=10))
                            env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end
                    end
                end
            end
            
            # Create data object for knapsack cuts
            dim_x = problem.n_facilities
            c_x = problem.fixed_costs
            dim_t = problem.n_customers # knapsack cut
            c_t = ones(dim_t)
            
            data = Data(dim_x, dim_t, problem, c_x, c_t)
            @assert dim_x == length(data.c_x)
            @assert dim_t == length(data.c_t)
            
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
                    lift = [true, false],
                    p in [1.0, Inf], 
                    disjunctive_cut_append_rule in [NoDisjunctiveCuts(), AllDisjunctiveCuts(), DisjunctiveCutsSmallerIndices()],
                    adjust_t_to_fx in [true; false]
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; lift $lift; p $p; dcut_append $disjunctive_cut_append_rule; adjust_t_to_fx $adjust_t_to_fx" begin
                        @info "solving UFLP p$i - disjunctive oracle/fatKnapsack - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"
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
                            lift = lift,
                            adjust_t_to_fx = adjust_t_to_fx
                        )
                        set_parameter!(disjunctive_oracle, oracle_param)
                        update_model!(disjunctive_oracle, data)
                        
                        @testset "NoSeq" begin
                            @info "solving UFLP p$i - disjunctive oracle/fatKnapsack/no seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"
                            master = Master(data; solver_param = master_solver_param)
                            update_model!(master, data)
                            root_preprocessing = NoRootNodePreprocessing()
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=10))
                            env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end
                        @testset "Seq" begin
                            @info "solving UFLP p$i - disjunctive oracle/fatKnapsack/seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"
                            master = Master(data; solver_param = master_solver_param)
                            update_model!(master, data)
                            root_preprocessing = RootNodePreprocessing(lazy_oracle, BendersSeq, BendersSeqParam(;time_limit=200.0, gap_tolerance=1e-6, verbose=false))
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=10))
                            env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end
                        @testset "SeqInOut" begin
                            @info "solving UFLP p$i - disjunctive oracle/fatKnapsack/seqinout - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"
                            master = Master(data; solver_param = master_solver_param)
                            update_model!(master, data)
                            root_preprocessing = RootNodePreprocessing(lazy_oracle, BendersSeqInOut, BendersSeqInOutParam(time_limit = 300.0, gap_tolerance = 1e-6, stabilizing_x = ones(data.dim_x), α = 0.9, λ = 0.1, verbose = false))
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=10))
                            env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end
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