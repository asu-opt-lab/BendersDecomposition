using BendersDecomposition
using Test
using JuMP
using CPLEX

@testset verbose = true "CFLP Disjunctive Callback Benders Tests" begin
    instances = setdiff(1:71, [67])

    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data
            data = read_cflp_benchmark_data("p$i")
            
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
            
            user_cb_param = UserCallbackParam(frequency=1)
            
            # Solve MIP for reference
            mip_model = Model()
            customize_mip_model!(mip_model, data)
            optimize!(mip_model)
            @assert termination_status(mip_model) == OPTIMAL
            mip_opt_val = objective_value(mip_model)
            
            @testset "Classic oracle" begin  
                # for strengthened in [true; false], add_benders_cuts_to_master in [true; false; 2], reuse_dcglp in [true; false], p in [1.0; Inf], lift in [true; false], disjunctive_cut_append_rule in [NoDisjunctiveCuts(); AllDisjunctiveCuts(); DisjunctiveCutsSmallerIndices()], adjust_t_to_fx in [true; false]
                for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [1.0], lift in [true], disjunctive_cut_append_rule in [AllDisjunctiveCuts()]
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; lift $lift; p $p; lift $lift; dcut_append $disjunctive_cut_append_rule" begin
                        oracle_param = SplitOracleParam(dcglp_param;
                                                                norm = LpNorm(p), 
                                                                split_index_selection_rule = RandomFractional(),
                                                                disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
                                                                strengthened = strengthened, 
                                                                add_benders_cuts_to_master = add_benders_cuts_to_master, 
                                                                fraction_of_benders_cuts_to_master = 0.5, 
                                                                reuse_dcglp = reuse_dcglp,
                                                                lift = lift)

                        @testset "NoSeq" begin
                            @info "solving CFLP p$i - disjunctive oracle/classical/no seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule"
                            master = Master(data; customize = customize_master_model!)
                            lazy_oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
                            typical_oracles = [ClassicalOracle(data, master; customize = customize_sub_model!), ClassicalOracle(data, master; customize = customize_sub_model!)]
                            disjunctive_oracle = SplitOracle(master, typical_oracles, oracle_param) 

                            root_preprocessing = NoRootNodePreprocessing()
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                            
                            env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                            if benders_param.verbose
                                @info "Disjunctive cuts added: $(length(env.user_callback.oracle.disjunctiveCuts))"
                                env.user_callback.oracle.param.add_benders_cuts_to_master != 0 && @info "Byproduct Benders cuts added: $(log.n_user_cuts[1] - length(env.user_callback.oracle.disjunctiveCuts))"
                            end
                        end

                        @testset "Seq" begin
                            @info "solving CFLP p$i - disjunctive oracle/classical/seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule"
                            master = Master(data; customize = customize_master_model!)
                            lazy_oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
                            typical_oracles = [ClassicalOracle(data, master; customize = customize_sub_model!), ClassicalOracle(data, master; customize = customize_sub_model!)]
                            disjunctive_oracle = SplitOracle(master, typical_oracles, oracle_param) 

                            root_preprocessing = RootNodePreprocessing(lazy_oracle, BendersSeq, BendersSeqParam(;time_limit=200.0, gap_tolerance=1e-9, verbose=false))
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                            
                            env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                            if benders_param.verbose
                                @info "Disjunctive cuts added: $(length(env.user_callback.oracle.disjunctiveCuts))"
                                env.user_callback.oracle.param.add_benders_cuts_to_master != 0 && @info "Byproduct Benders cuts added: $(log.n_user_cuts[1] - length(env.user_callback.oracle.disjunctiveCuts))"
                            end
                        end

                        @testset "SeqInOut" begin
                            @info "solving CFLP p$i - disjunctive oracle/classical/seqinout - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule"
                            master = Master(data; customize = customize_master_model!)
                            lazy_oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
                            typical_oracles = [ClassicalOracle(data, master; customize = customize_sub_model!), ClassicalOracle(data, master; customize = customize_sub_model!)]
                            disjunctive_oracle = SplitOracle(master, typical_oracles, oracle_param) 

                            root_preprocessing = RootNodePreprocessing(lazy_oracle, BendersSeqInOut, BendersSeqInOutParam(time_limit = 300.0, gap_tolerance = 1e-9, stabilizing_x = ones(data.n_facilities), α = 0.9, λ = 0.1, verbose = false))
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                            
                            env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                            if benders_param.verbose
                                @info "Disjunctive cuts added: $(length(env.user_callback.oracle.disjunctiveCuts))"
                                env.user_callback.oracle.param.add_benders_cuts_to_master != 0 && @info "Byproduct Benders cuts added: $(log.n_user_cuts[1] - length(env.user_callback.oracle.disjunctiveCuts))"
                            end
                        end
                    end
                end
            end
            
            @testset "Knapsack oracle" begin
                # for strengthened in [true; false], add_benders_cuts_to_master in [true; false; 2], reuse_dcglp in [true; false], p in [1.0; Inf], lift in [true; false], disjunctive_cut_append_rule in [NoDisjunctiveCuts(); AllDisjunctiveCuts(); DisjunctiveCutsSmallerIndices()], adjust_t_to_fx in [true; false]
                for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [1.0], lift in [true], disjunctive_cut_append_rule in [AllDisjunctiveCuts()]
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; lift $lift; p $p; dcut_append $disjunctive_cut_append_rule" begin
                        
                        oracle_param = SplitOracleParam(dcglp_param;
                                                                norm = LpNorm(p),
                                                                split_index_selection_rule = RandomFractional(),
                                                                disjunctive_cut_append_rule = disjunctive_cut_append_rule,
                                                                strengthened = strengthened,
                                                                add_benders_cuts_to_master = add_benders_cuts_to_master,
                                                                fraction_of_benders_cuts_to_master = 0.5,
                                                                reuse_dcglp = reuse_dcglp,
                                                                lift = lift)

                        @testset "NoSeq" begin
                            @info "solving CFLP p$i - disjunctive oracle/knapsack oracle/no seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule"
                            master = Master(data; customize = customize_master_model!)
                            lazy_oracle = CFLKnapsackOracle(data, master; customize = customize_sub_model!)
                            typical_oracles = [CFLKnapsackOracle(data, master; customize = customize_sub_model!), CFLKnapsackOracle(data, master; customize = customize_sub_model!)]
                            disjunctive_oracle = SplitOracle(master, typical_oracles, oracle_param) 

                            root_preprocessing = NoRootNodePreprocessing()
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                            
                            env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                            if benders_param.verbose
                                @info "Disjunctive cuts added: $(length(env.user_callback.oracle.disjunctiveCuts))"
                                env.user_callback.oracle.param.add_benders_cuts_to_master != 0 && @info "Byproduct Benders cuts added: $(log.n_user_cuts[1] - length(env.user_callback.oracle.disjunctiveCuts))"
                            end
                        end
                        @testset "Seq" begin
                            @info "solving CFLP p$i - disjunctive oracle/knapsack oracle/seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule"
                            master = Master(data; customize = customize_master_model!)
                            lazy_oracle = CFLKnapsackOracle(data, master; customize = customize_sub_model!)
                            typical_oracles = [CFLKnapsackOracle(data, master; customize = customize_sub_model!), CFLKnapsackOracle(data, master; customize = customize_sub_model!)]
                            disjunctive_oracle = SplitOracle(master, typical_oracles, oracle_param) 
                            
                            root_preprocessing = RootNodePreprocessing(lazy_oracle, BendersSeq, BendersSeqParam(;time_limit=200.0, gap_tolerance=1e-9, verbose=false))
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                            
                            env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                            if benders_param.verbose
                                @info "Disjunctive cuts added: $(length(env.user_callback.oracle.disjunctiveCuts))"
                                env.user_callback.oracle.param.add_benders_cuts_to_master != 0 && @info "Byproduct Benders cuts added: $(log.n_user_cuts[1] - length(env.user_callback.oracle.disjunctiveCuts))"
                            end
                        end
                        @testset "SeqInOut" begin
                            @info "solving CFLP p$i - disjunctive oracle/knapsack oracle/seqinout - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule"
                            master = Master(data; customize = customize_master_model!)
                            lazy_oracle = CFLKnapsackOracle(data, master; customize = customize_sub_model!)
                            typical_oracles = [CFLKnapsackOracle(data, master; customize = customize_sub_model!), CFLKnapsackOracle(data, master; customize = customize_sub_model!)]
                            disjunctive_oracle = SplitOracle(master, typical_oracles, oracle_param) 

                            root_preprocessing = RootNodePreprocessing(lazy_oracle, BendersSeqInOut, BendersSeqInOutParam(time_limit = 300.0, gap_tolerance = 1e-9, stabilizing_x = ones(data.n_facilities), α = 0.9, λ = 0.1, verbose = false))
                            lazy_callback = LazyCallback(lazy_oracle)
                            user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                            
                            env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                            if benders_param.verbose
                                @info "Disjunctive cuts added: $(length(env.user_callback.oracle.disjunctiveCuts))"
                                env.user_callback.oracle.param.add_benders_cuts_to_master != 0 && @info "Byproduct Benders cuts added: $(log.n_user_cuts[1] - length(env.user_callback.oracle.disjunctiveCuts))"
                            end
                        end
                    end
                end
            end
        end
    end
end