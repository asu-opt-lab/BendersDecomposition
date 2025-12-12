# tight tolerance for dcglp; adjust tol for generate_cut with omega_0
# dcglp is not solved -> cannot proceed the algorithm, and should be terminated. Added an optimal argument to dcglp `throw_typical_cuts_for_errors=false`
# remove other_constraints from classical oracle.
# add a method for generating optimal vertex for SpecializedBendersSeq

using BendersDecomposition
using Test
using JuMP
using CPLEX
using Logging
global_logger(ConsoleLogger(stderr, Logging.Debug))

# loop parameters
specialized_benders_param = SpecializedBendersSeqParam(;
                                                        time_limit = 1000.0,
                                                        lp_gap_tolerance = 1e-9,
                                                        integrality_tolerance = 1e-9,
                                                        verbose = true
                                                        )
dcglp_optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_LPMETHOD" => 3, MOI.Silent() => true)
dcglp_param = DcglpParam(dcglp_optimizer;
                        time_limit = 1000.0, 
                        gap_tolerance = 1e-9, 
                        halt_limit = Int(1e9), 
                        iter_limit = Int(1e9),
                        verbose = true
                        )

@testset verbose = true "UFLP Specialized Sequential Benders Tests" begin
    # instances = setdiff(1:71, [67])
    instances = [49] # only instance that uses SpecializedBendersSeq
    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data if necessary
            problem = read_uflp_benchmark_data("p$(i)")
            
            # Solve MIP for reference
            mip_model = Model()
            customize_mip_model!(mip_model, problem)
            optimize!(mip_model)
            @assert termination_status(mip_model) == OPTIMAL
            mip_opt_val = objective_value(mip_model)

            @testset "Classical oracle" begin
                # for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf]
                # for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [Inf] # too slow

                for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0]
                # for strengthened in [false], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0]
                # for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [1.0]
                # for strengthened in [false], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [1.0]

                # for strengthened in [false], add_benders_cuts_to_master in [false], reuse_dcglp in [false], p in [1.0]
                # for strengthened in [true], add_benders_cuts_to_master in [false], reuse_dcglp in [false], p in [1.0]
                # for strengthened in [true], add_benders_cuts_to_master in [false], reuse_dcglp in [true], p in [1.0] # fail
                # for strengthened in [false], add_benders_cuts_to_master in [false], reuse_dcglp in [true], p in [1.0] #fail
                    @info "solving p$i - begin oracle - Specialized seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p"
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p" begin
                        oracle_param = DisjunctiveOracleParam(dcglp_param;
                                                                norm = LpNorm(p), 
                                                                split_index_selection_rule = LargestFractional(),
                                                                disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(), 
                                                                strengthened=strengthened, 
                                                                add_benders_cuts_to_master=add_benders_cuts_to_master, 
                                                                fraction_of_benders_cuts_to_master = 0.5, 
                                                                reuse_dcglp=reuse_dcglp)

                        master = Master(problem; customize = customize_master_model!)
                        set_optimizer_attribute(master.model, "CPX_PARAM_LPMETHOD", 1)
                        typical_oracles = [ClassicalOracle(problem, master; customize = customize_sub_model!); ClassicalOracle(problem, master; customize = customize_sub_model!)] # for kappa & nu
                        disjunctive_oracle = DisjunctiveOracle(master, typical_oracles, oracle_param) 
                        env = SpecializedBendersSeq(master, disjunctive_oracle; param = specialized_benders_param)
                        
                        log = solve!(env)
                        @test env.termination_status == Optimal() ? isapprox(mip_opt_val, env.obj_value, atol=1e-5) : false
                    end
                end
            end 

            @testset "fat knapsack oracle" begin
                function customize_master_model!(model::Model, problem::UFLPData)
                    optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)
                    set_optimizer(model, optimizer)
                    @variable(model, x[1:problem.n_facilities], Bin)
                    @variable(model, t[1:problem.n_customers] >= -1e6)
                    @constraint(model, sum(x) >= 2)
                    @objective(model, Min, problem.fixed_costs'* x + sum(t))
                    return (x = x, ), t
                end

                for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0] 
                # for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf] 
                    @info "solving p$i - fat Knapsack oracle - Specialized seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p"
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p" begin
                        oracle_param = DisjunctiveOracleParam(dcglp_param;
                                                                norm = LpNorm(p), 
                                                                split_index_selection_rule = LargestFractional(),
                                                                disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(), 
                                                                strengthened=strengthened, 
                                                                add_benders_cuts_to_master=add_benders_cuts_to_master, 
                                                                fraction_of_benders_cuts_to_master = 0.5, 
                                                                reuse_dcglp=reuse_dcglp)

                        master = Master(problem; customize = customize_master_model!)
                        set_optimizer_attribute(master.model, "CPX_PARAM_LPMETHOD", 1)
                        typical_oracles = [UFLKnapsackOracle(problem); UFLKnapsackOracle(problem)] # for kappa & nu
                        disjunctive_oracle = DisjunctiveOracle(master, typical_oracles, oracle_param) 
                        env = SpecializedBendersSeq(master, disjunctive_oracle; param = specialized_benders_param)
                        
                        log = solve!(env)
                        @test env.termination_status == Optimal() ? isapprox(mip_opt_val, env.obj_value, atol=1e-5) : false
                    end
                end
            end
        end
    end
end

@testset verbose = true "CFLP Specialized Sequential Benders Tests" begin
    # instances = setdiff(1:71, [67])
    instances = [25 32 34 36 49 51]
    # numerical issue: 29 30 31 35
    # success: 25 32 34 36 49 51
    # lp: 1:24, 26:28, 33, 37-48, 50, 52-66, 68-71
    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data if necessary
            problem = read_cflp_benchmark_data("p$i")
            # problem = read_GK_data("f100-c100-r3-1")
                            
            # Solve MIP for reference
            mip_model = Model()
            customize_mip_model!(mip_model, problem)
            optimize!(mip_model)
            @assert termination_status(mip_model) == OPTIMAL
            mip_opt_val = objective_value(mip_model)
            x_opt = value.(mip_model[:x])
            t_opt = value.(mip_model[:t])
            @info "MIP_opt_val: $mip_opt_val"
            @info "MIP sol: $(x_opt)"

            @testset "Classical oracle" begin 
                    # for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf]
                for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0]
                    @info "solving CFLP p$i - disjunctive oracle/classical - specialized seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p"
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p" begin
                        oracle_param = DisjunctiveOracleParam(dcglp_param;
                                                                norm = LpNorm(p), 
                                                                split_index_selection_rule = LargestFractional(),
                                                                disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(), 
                                                                strengthened=strengthened, 
                                                                add_benders_cuts_to_master=add_benders_cuts_to_master, 
                                                                fraction_of_benders_cuts_to_master = 0.5, 
                                                                reuse_dcglp=reuse_dcglp)

                        master = Master(problem; customize = customize_master_model!)
                        typical_oracles = [ClassicalOracle(problem, master; customize = customize_sub_model!); ClassicalOracle(problem, master; customize = customize_sub_model!)] # for kappa & nu
                        disjunctive_oracle = DisjunctiveOracle(master, typical_oracles, oracle_param) 
                        env = SpecializedBendersSeq(master, disjunctive_oracle; param = specialized_benders_param)
                        
                        log = solve!(env)
                        @test env.termination_status == Optimal()
                        # if env.log.termination_status == Optimal()
                        if !isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                            t_opt_ = [sum(t_opt)]
                            @error "Failed ***** mip_opt_val = $(mip_opt_val) vs BD_obj_val = $(env.obj_value)"
                            # optimize!(master.model)
                            opt_sol = Dict{VariableRef, Float64}()
                            for i = 1:master.dim_x
                                opt_sol[master.x[i]] = x_opt[i]
                            end
                            for i = 1:master.dim_t
                                opt_sol[master.t[i]] = t_opt_[i]
                            end

                            @info primal_feasibility_report(env.master.model, opt_sol)
                            @info master.c_x' * x_opt + master.c_t' * t_opt_
                            
                            for v in keys(opt_sol)
                                fix(v, opt_sol[v]; force=true)
                            end
                            optimize!(env.master.model)
                            @info objective_value(env.master.model)
                        end
                        @test env.termination_status == Optimal() ? isapprox(mip_opt_val, env.obj_value, atol=1e-5) : false
                    end
                end
            end 

            @testset "Knapsack oracle" begin
                # for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf] 
                for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0] 
                    @info "solving CFLP p$i - disjunctive oracle/knapsack- specialized seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p"
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p" begin
                            oracle_param = DisjunctiveOracleParam(dcglp_param;
                                                                norm = LpNorm(p), 
                                                                split_index_selection_rule = LargestFractional(),
                                                                disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(), 
                                                                strengthened=strengthened, 
                                                                add_benders_cuts_to_master=add_benders_cuts_to_master, 
                                                                fraction_of_benders_cuts_to_master = 0.5, 
                                                                reuse_dcglp=reuse_dcglp)
                            
                                                                master = Master(problem; customize = customize_master_model!)
                            typical_oracles = [CFLKnapsackOracle(problem, master; customize = customize_sub_model!); CFLKnapsackOracle(problem, master; customize = customize_sub_model!)]
                            disjunctive_oracle = DisjunctiveOracle(master, typical_oracles, oracle_param) 
                            env = SpecializedBendersSeq(master, disjunctive_oracle; param = specialized_benders_param)
                        log = solve!(env)
                        @test env.termination_status == Optimal() ? isapprox(mip_opt_val, env.obj_value, atol=1e-5) : false
                    end
                end
            end
        end
    end
end