include("$(dirname(dirname(@__DIR__)))/example/snip/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/snip/model.jl")

@testset verbose = true "SNIP Sequential Benders Tests" begin
    for instance in [0], snipno in [0], budget in [30.0]
        @testset "instance $instance; snipno $snipno budget $budget" begin
            # Load problem data
            problem = read_snip_data(instance, snipno, budget)
            
            # Initialize data object
            dim_x = length(problem.D)
            dim_t = problem.num_scenarios
            c_x = zeros(dim_x)
            c_t = map(p -> p[3], problem.scenarios)
            data = Data(dim_x, dim_t, problem, c_x, c_t)
            @assert dim_x == length(data.c_x)
            @assert dim_t == length(data.c_t)
            
            # Get standard parameters
            benders_param = BendersBnBParam(;
                time_limit = 3600.0,
                gap_tolerance = 1e-6,
                verbose = true
            )

            dcglp_param = DcglpParam(
                time_limit = 200.0,
                gap_tolerance = 1e-3,
                halt_limit = 3,
                iter_limit = 15,
                verbose = true
            )

            user_cb_param = UserCallbackParam(frequency=100)
            
            # Common solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPXPARAM_Threads" => 4)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 4)
            typical_oracle_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)
            dcglp_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9) 
            
            # Solve MIP for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            set_optimizer_attribute(mip.model, MOI.Silent(), false)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)
            set_optimizer_attribute(mip.model, MOI.Silent(), true)

            mip_opt_val = 0.3131350534615714
            
            @testset "Classic oracle" begin
                @testset "NoSeq" begin

                    lazy_oracle = SeparableOracle(data, ClassicalOracle(), data.problem.num_scenarios; solver_param = typical_oracle_solver_param)
                    for j=1:lazy_oracle.N
                        update_model!(lazy_oracle.oracles[j], data, j)
                    end
                    typical_oracle_kappa = SeparableOracle(data, ClassicalOracle(), data.problem.num_scenarios; solver_param = typical_oracle_solver_param)
                    typical_oracle_nu = SeparableOracle(data, ClassicalOracle(), data.problem.num_scenarios; solver_param = typical_oracle_solver_param)
                    for j=1:typical_oracle_kappa.N
                        update_model!(typical_oracle_kappa.oracles[j], data, j)
                    end
                    for j=1:typical_oracle_nu.N
                        update_model!(typical_oracle_nu.oracles[j], data, j)
                    end
                    typical_oracles = [typical_oracle_kappa; typical_oracle_nu]
                
                    # Test various parameter combinations
                    for strengthened in [true, false], 
                        add_benders_cuts_to_master in [true; 2], 
                        reuse_dcglp in [true, false], 
                        lift in [true, false],
                        p in [1.0, Inf], 
                        disjunctive_cut_append_rule in [NoDisjunctiveCuts(), AllDisjunctiveCuts(), DisjunctiveCutsSmallerIndices()],
                        adjust_t_to_fx in [true; false]
                        
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; lift $lift; p $p; dcut_append $disjunctive_cut_append_rule;  adjust_t_to_fx $adjust_t_to_fx" begin
                            disjunctive_oracle = DisjunctiveOracle(data, typical_oracles; solver_param = dcglp_solver_param, param = dcglp_param) 

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
                                @info "solving SNIP instance$instance snipno $snipno budget $budget - disjunctive oracle/classical/no seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"
                                master = Master(data; solver_param = master_solver_param)
                                update_model!(master, data)
                                root_preprocessing = NoRootNodePreprocessing()
                                lazy_callback = LazyCallback(lazy_oracle)
                                user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                                env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                                log = solve!(env)
                                @test env.termination_status == Optimal()
                                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                            end

                            # @testset "Seq" begin
                            #     @info "solving SNIP instance$instance snipno $snipno budget $budget - disjunctive oracle/classical/seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"
                            #     master = Master(data; solver_param = master_solver_param)
                            #     update_model!(master, data)
                            #     root_preprocessing = RootNodePreprocessing(lazy_oracle, BendersSeq, BendersSeqParam(;time_limit=200.0, gap_tolerance=1e-6, verbose=true))
                            #     lazy_callback = LazyCallback(lazy_oracle)
                            #     user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                            #     env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            #     log = solve!(env)
                            #     @test env.termination_status == Optimal()
                            #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                            # end

                            # @testset "SeqInOut" begin
                            #     @info "solving SNIP instance$instance snipno $snipno budget $budget - disjunctive oracle/classical/seqinout - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp lift $lift p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"
                            #     master = Master(data; solver_param = master_solver_param)
                            #     update_model!(master, data)
                            #     root_preprocessing = RootNodePreprocessing(lazy_oracle, BendersSeqInOut, BendersSeqInOutParam(time_limit = 300.0, gap_tolerance = 1e-6, stabilizing_x = ones(data.dim_x), α = 0.9, λ = 0.1, verbose = false))
                            #     lazy_callback = LazyCallback(lazy_oracle)
                            #     user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
                            #     env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                            #     log = solve!(env)
                            #     @test env.termination_status == Optimal()
                            #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                            # end
                        end
                    end
                end
            end

        end
    end
end


            
            
            
            
            



