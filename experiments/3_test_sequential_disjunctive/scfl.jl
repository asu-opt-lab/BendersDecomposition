using BendersDecomposition
using Test
using JuMP
using CPLEX

@testset verbose = true "Stochastic CFLP Sequential Benders Tests" begin
    instances = 1:5

    for i in instances
        @testset "Instance: f25-c50-s64-r10-$i" begin
            # Load problem data
            problem = read_stochastic_capacited_facility_location_problem("f25-c50-s64-r10-$i")
            
            # Algorithm parameters
            benders_param = BendersSeqParam(;
                                            time_limit = 2000.0,
                                            gap_tolerance = 1e-6,
                                            verbose = false
                                            )
            dcglp_optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9, MOI.Silent() => true)
            dcglp_param = DcglpParam(dcglp_optimizer;
                                    time_limit = 1000.0, 
                                    gap_tolerance = 1e-3, 
                                    halt_limit = 3, 
                                    iter_limit = 250,
                                    verbose = false
                                    )
            
            # Solve MIP for reference
            mip_model = Model()
            customize_mip_model!(mip_model, problem)
            optimize!(mip_model)
            @assert termination_status(mip_model) == OPTIMAL
            mip_opt_val = objective_value(mip_model)

            @testset "Classic oracle" begin
                @testset "Seq" begin        
                    # for strengthened in [true; false], add_benders_cuts_to_master in [true; false; 2], reuse_dcglp in [true; false], p in [1.0; Inf], lift in [true; false], disjunctive_cut_append_rule in [NoDisjunctiveCuts(); AllDisjunctiveCuts(); DisjunctiveCutsSmallerIndices()], adjust_t_to_fx in [true; false]
                    for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [1.0], lift in [true], disjunctive_cut_append_rule in [AllDisjunctiveCuts()]
                        @info "solving SCFLP f25-c50-s64-r10-$i - disjunctive oracle/classical - seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p lift $lift dcut_append $disjunctive_cut_append_rule"
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; p $p; lift $lift; dcut_append $disjunctive_cut_append_rule" begin
                            
                            oracle_param = DisjunctiveOracleParam(dcglp_param;
                                                                norm = LpNorm(p), 
                                                                split_index_selection_rule = RandomFractional(),
                                                                disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
                                                                strengthened = strengthened, 
                                                                add_benders_cuts_to_master = add_benders_cuts_to_master, 
                                                                fraction_of_benders_cuts_to_master = 1.0, 
                                                                reuse_dcglp = reuse_dcglp,
                                                                lift = lift)
                                                                    
                            master = Master(problem; customize = customize_master_model!)
                            typical_oracle_kappa = SeparableOracle(problem, master, ClassicalOracle(), problem.n_scenarios; customize = customize_sub_model!)
                            typical_oracle_nu = SeparableOracle(problem, master, ClassicalOracle(), problem.n_scenarios; customize = customize_sub_model!)
                            typical_oracles = [typical_oracle_kappa; typical_oracle_nu]
                            disjunctive_oracle = DisjunctiveOracle(master, typical_oracles, oracle_param) 
                            env = BendersSeq(master, disjunctive_oracle; param = benders_param)

                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            # if env.log.termination_status == Optimal()
                            if !isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                                t_opt_ = [sum(t_opt)]
                                @error "Failed ***** mip_opt_val = $(mip_opt_val) vs BD_obj_val = $(env.obj_value)" 
                                # optimize!(master.model)
                                opt_sol = Dict{VariableRef, Float64}()
                                for i = 1:data.dim_x
                                    opt_sol[master.x[i]] = x_opt[i]
                                end
                                for i = 1:data.dim_t
                                    opt_sol[master.t[i]] = t_opt_[i]
                                end

                                @info primal_feasibility_report(env.master.model, opt_sol)
                                @info data.c_x' * x_opt + data.c_t' * t_opt_

                                for v in keys(opt_sol)
                                    fix(v, opt_sol[v]; force=true)
                                end
                                optimize!(env.master.model)
                                @info objective_value(env.master.model)
                            end
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                            # elseif env.termination_status == TimeLimit()
                            #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                            #     @test env.log.LB <= mip_opt_val <= env.log.UB
                            # elseif env.log.termination_status == InfeasibleOrNumericalIssue()
                            #     @test false
                            # end
                        end
                    end
                end
            end 
            
            @testset "Knapsack oracle" begin
                @testset "Seq" begin        
                    # for strengthened in [true; false], add_benders_cuts_to_master in [true; false; 2], reuse_dcglp in [true; false], p in [1.0; Inf], lift in [true; false], disjunctive_cut_append_rule in [NoDisjunctiveCuts(); AllDisjunctiveCuts(); DisjunctiveCutsSmallerIndices()], adjust_t_to_fx in [true; false]
                    for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [1.0], lift in [true], disjunctive_cut_append_rule in [AllDisjunctiveCuts()]
                        @info "solving SCFLP f25-c50-s64-r10-$i - disjunctive oracle/classical - seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p lift $lift dcut_append $disjunctive_cut_append_rule"
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master; reuse $reuse_dcglp; p $p; lift $lift; dcut_append $disjunctive_cut_append_rule" begin

                            oracle_param = DisjunctiveOracleParam(dcglp_param;
                                                                norm = LpNorm(p), 
                                                                split_index_selection_rule = RandomFractional(),
                                                                disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
                                                                strengthened = strengthened, 
                                                                add_benders_cuts_to_master = add_benders_cuts_to_master, 
                                                                fraction_of_benders_cuts_to_master = 1.0, 
                                                                reuse_dcglp = reuse_dcglp,
                                                                lift = lift)
                                                                    
                            master = Master(problem; customize = customize_master_model!)
                            typical_oracle_kappa = SeparableOracle(problem, master, CFLKnapsackOracle(), problem.n_scenarios; customize = customize_sub_model!)
                            typical_oracle_nu = SeparableOracle(problem, master, CFLKnapsackOracle(), problem.n_scenarios; customize = customize_sub_model!)
                            typical_oracles = [typical_oracle_kappa; typical_oracle_nu]
                            disjunctive_oracle = DisjunctiveOracle(master, typical_oracles, oracle_param) 
                            env = BendersSeq(master, disjunctive_oracle; param = benders_param)

                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            # if env.log.termination_status == Optimal()
                            if !isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                                t_opt_ = [sum(t_opt)]
                                @error "Failed ***** mip_opt_val = $(mip_opt_val) vs BD_obj_val = $(env.obj_value)" 
                                # optimize!(master.model)
                                opt_sol = Dict{VariableRef, Float64}()
                                for i = 1:data.dim_x
                                    opt_sol[master.x[i]] = x_opt[i]
                                end
                                for i = 1:data.dim_t
                                    opt_sol[master.t[i]] = t_opt_[i]
                                end

                                @info primal_feasibility_report(env.master.model, opt_sol)
                                @info data.c_x' * x_opt + data.c_t' * t_opt_

                                for v in keys(opt_sol)
                                    fix(v, opt_sol[v]; force=true)
                                end
                                optimize!(env.master.model)
                                @info objective_value(env.master.model)
                            end
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                            # elseif env.termination_status == TimeLimit()
                            #     @warn "TIME LIMIT, elapsed time = $(time() - env.log.start_time)"
                            #     @test env.log.LB <= mip_opt_val <= env.log.UB
                            # elseif env.log.termination_status == InfeasibleOrNumericalIssue()
                            #     @test false
                            # end
                        end
                    end
                end
            end
        end
    end
end