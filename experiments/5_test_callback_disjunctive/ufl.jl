using BendersDecomposition
using Test
using JuMP
using CPLEX

@testset verbose = true "UFLP Callback Disjunctive Benders Tests" begin
    instances = setdiff(1:71, [67])

    for i in instances
        @testset "Instance: p$i" begin            
            # Load problem data
            data = read_uflp_benchmark_data("p$i")

            # Algorithm parameters
            benders_param = BendersBnBParam(;
                                            time_limit = 200.0,
                                            gap_tolerance = 1e-6,
                                            verbose = false
                                            )

            dcglp_optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9, MOI.Silent() => true)
            dcglp_param = DcglpParam(dcglp_optimizer;
                                    time_limit = 200.0, 
                                    gap_tolerance = 1e-3, 
                                    halt_limit = 3, 
                                    iter_limit = 15,
                                    verbose = false
                                    )
            
            # Solve MIP for reference
            mip_model = Model()
            customize_mip_model!(mip_model, data)
            set_optimizer_attribute(mip_model, "CPX_PARAM_BRDIR", 1)
            optimize!(mip_model)
            @assert termination_status(mip_model) == OPTIMAL
            mip_opt_val = objective_value(mip_model)
            
            # Test classical oracle
            @testset "Classic oracle" begin
                # for strengthened in [true; false], add_benders_cuts_to_master in [true; false; 2], reuse_dcglp in [true; false], p in [1.0; Inf], lift in [true; false], disjunctive_cut_append_rule in [NoDisjunctiveCuts(); AllDisjunctiveCuts(); DisjunctiveCutsSmallerIndices()], adjust_t_to_fx in [true; false]
                for strengthened in [true], add_benders_cuts_to_master in [2], reuse_dcglp in [false], p in [1.0], lift in [false], disjunctive_cut_append_rule in [AllDisjunctiveCuts()]
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; p $p; lift $lift; dcut_append $disjunctive_cut_append_rule" begin
                        @info "solving UFLP p$i - disjunctive oracle/classical - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p lift $lift dcut_append $disjunctive_cut_append_rule"   
                        
                        oracle_param = SplitOracleParam(dcglp_param;
                                                            norm = LpNorm(p), 
                                                            split_index_selection_rule = LargestFractional(),
                                                            disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
                                                            strengthened = strengthened, 
                                                            add_benders_cuts_to_master = add_benders_cuts_to_master, 
                                                            fraction_of_benders_cuts_to_master = 0.05, 
                                                            reuse_dcglp = reuse_dcglp,
                                                            lift = lift)

                        @testset "NoSeq" begin
                            @info "solving UFLP p$i - disjunctive oracle/classical/no seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule"
                            master = Master(data; customize = customize_master_model!)
                            set_optimizer_attribute(master.model, "CPX_PARAM_BRDIR", 1)
                            lazy_oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
                            typical_oracles = [ClassicalOracle(data, master; customize = customize_sub_model!); ClassicalOracle(data, master; customize = customize_sub_model!)] # for kappa & nu
                            disjunctive_oracle = SplitOracle(master, typical_oracles, oracle_param) 
                    
                            root_preprocessing = NoRootNodePreprocessing()
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=10))
                            
                            env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end

                        @testset "Seq" begin
                            @info "solving UFLP p$i - disjunctive oracle/classical/seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule"
                            master = Master(data; customize = customize_master_model!)
                            set_optimizer_attribute(master.model, "CPX_PARAM_BRDIR", 1)
                            lazy_oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
                            typical_oracles = [ClassicalOracle(data, master; customize = customize_sub_model!); ClassicalOracle(data, master; customize = customize_sub_model!)] # for kappa & nu
                            disjunctive_oracle = SplitOracle(master, typical_oracles, oracle_param) 
                    
                            root_preprocessing = RootNodePreprocessing(lazy_oracle, BendersSeq, BendersSeqParam(;time_limit=200.0, gap_tolerance=1e-9, verbose=false))
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=10))
                            
                            env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end

                        @testset "SeqInOut" begin
                            @info "solving UFLP p$i - disjunctive oracle/classical/seqinout - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule"
                            master = Master(data; customize = customize_master_model!)
                            set_optimizer_attribute(master.model, "CPX_PARAM_BRDIR", 1)
                            lazy_oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
                            typical_oracles = [ClassicalOracle(data, master; customize = customize_sub_model!); ClassicalOracle(data, master; customize = customize_sub_model!)] # for kappa & nu
                            disjunctive_oracle = SplitOracle(master, typical_oracles, oracle_param) 

                            root_preprocessing = RootNodePreprocessing(lazy_oracle, BendersSeqInOut, BendersSeqInOutParam(time_limit = 300.0, gap_tolerance = 1e-9, stabilizing_x = ones(data.n_facilities), α = 0.9, λ = 0.1, verbose = false))
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=10))
                            
                            env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end
                    end
                end
            end
            
            @testset "Fat knapsack oracle" begin
                function customize_master_model!(model::Model, data::UFLPData)
                    optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)
                    set_optimizer(model, optimizer)
                    @variable(model, x[1:data.n_facilities], Bin)
                    @variable(model, t[1:data.n_customers] >= -1e6)
                    @constraint(model, sum(x) >= 2)
                    @objective(model, Min, data.fixed_costs'* x + sum(t))
                    return (x = x, ), t
                end

                # for strengthened in [true; false], add_benders_cuts_to_master in [true; false; 2], reuse_dcglp in [true; false], p in [1.0; Inf], lift in [true; false], disjunctive_cut_append_rule in [NoDisjunctiveCuts(); AllDisjunctiveCuts(); DisjunctiveCutsSmallerIndices()], adjust_t_to_fx in [true; false]
                for strengthened in [true], add_benders_cuts_to_master in [2], reuse_dcglp in [false], p in [1.0], lift in [false], disjunctive_cut_append_rule in [AllDisjunctiveCuts()]
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; p $p; lift $lift; dcut_append $disjunctive_cut_append_rule" begin
                        @info "solving UFLP p$i - disjunctive oracle/fat knapsack - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p lift $lift dcut_append $disjunctive_cut_append_rule"
                        
                        oracle_param = SplitOracleParam(dcglp_param;
                                                            norm = LpNorm(p), 
                                                            split_index_selection_rule = LargestFractional(),
                                                            disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
                                                            strengthened = strengthened, 
                                                            add_benders_cuts_to_master = add_benders_cuts_to_master, 
                                                            fraction_of_benders_cuts_to_master = 0.05, 
                                                            reuse_dcglp = reuse_dcglp,
                                                            lift = lift)
                        
                        @testset "NoSeq" begin
                            @info "solving UFLP p$i - disjunctive oracle/fat knapsack/no seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule"
                            master = Master(data; customize = customize_master_model!)
                            set_optimizer_attribute(master.model, "CPX_PARAM_BRDIR", 1)
                            lazy_oracle = UFLKnapsackOracle(data)
                            set_parameter!(lazy_oracle, "add_only_violated_cuts", true)
                            typical_oracles = [UFLKnapsackOracle(data), UFLKnapsackOracle(data)]
                            for k=1:2
                                set_parameter!(typical_oracles[k], "add_only_violated_cuts", true)
                            end
                            disjunctive_oracle = SplitOracle(master, typical_oracles, oracle_param) 

                            root_preprocessing = NoRootNodePreprocessing()
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=10))
                            
                            env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end

                        @testset "Seq" begin
                            @info "solving UFLP p$i - disjunctive oracle/fat knapsack/seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule"
                            master = Master(data; customize = customize_master_model!)
                            set_optimizer_attribute(master.model, "CPX_PARAM_BRDIR", 1)
                            lazy_oracle = UFLKnapsackOracle(data)
                            set_parameter!(lazy_oracle, "add_only_violated_cuts", true)
                            typical_oracles = [UFLKnapsackOracle(data), UFLKnapsackOracle(data)]
                            for k=1:2
                                set_parameter!(typical_oracles[k], "add_only_violated_cuts", true)
                            end
                            disjunctive_oracle = SplitOracle(master, typical_oracles, oracle_param) 

                            root_preprocessing = RootNodePreprocessing(lazy_oracle, BendersSeq, BendersSeqParam(;time_limit=200.0, gap_tolerance=1e-9, verbose=false))
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=10))
                            
                            env = BendersBnB( master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end

                        @testset "SeqInOut" begin
                            @info "solving UFLP p$i - disjunctive oracle/fat knapsack/seqinout - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule"
                            master = Master(data; customize = customize_master_model!)
                            set_optimizer_attribute(master.model, "CPX_PARAM_BRDIR", 1)
                            lazy_oracle = UFLKnapsackOracle(data)
                            set_parameter!(lazy_oracle, "add_only_violated_cuts", true)
                            typical_oracles = [UFLKnapsackOracle(data), UFLKnapsackOracle(data)]
                            for k=1:2
                                set_parameter!(typical_oracles[k], "add_only_violated_cuts", true)
                            end
                            disjunctive_oracle = SplitOracle(master, typical_oracles, oracle_param) 

                            root_preprocessing = RootNodePreprocessing(lazy_oracle, BendersSeqInOut, BendersSeqInOutParam(time_limit = 300.0, gap_tolerance = 1e-9, stabilizing_x = ones(data.n_facilities), α = 0.9, λ = 0.1, verbose = false))
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=10))
                            
                            env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                        end
                    end
                end
            end
        end
    end
end