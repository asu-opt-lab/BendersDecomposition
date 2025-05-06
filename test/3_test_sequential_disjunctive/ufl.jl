include("$(dirname(dirname(@__DIR__)))/example/uflp/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/uflp/oracle.jl")
include("$(dirname(dirname(@__DIR__)))/example/uflp/model.jl")

@testset verbose = true "UFLP Sequential Benders Tests -- MIP master" begin
    instances = setdiff(1:71, [67])
    # instances = 49:49
    # instances = setdiff(41:71, [67])
    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data if necessary
            problem = read_uflp_benchmark_data("p$(i)")
            
            # initialize dim_x, dim_t, c_x, c_t
            dim_x = problem.n_facilities
            c_x = problem.fixed_costs
            dim_t = 1 # classical cut
            c_t = [1]
            
            data = Data(dim_x, dim_t, problem, c_x, c_t)
            @assert dim_x == length(data.c_x)
            @assert dim_t == length(data.c_t)

            # algorithm parameters
            benders_param = BendersSeqParam(;
                            time_limit = 800.0,
                            gap_tolerance = 1e-6,
                            verbose = false
                        )
            benders_inout_param = BendersSeqInOutParam(;
                                        time_limit = 800.0,
                                        gap_tolerance = 1e-6,
                                        verbose = false,
                                        stabilizing_x = ones(data.dim_x),
                                        α = 0.9,
                                        λ = 0.1
                                    )
            dcglp_param = DcglpParam(;
                                    time_limit = 1000.0, 
                                    gap_tolerance = 1e-3, 
                                    halt_limit = 3, 
                                    iter_limit = 250,
                                    verbose = false
                            )

            # solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 4)
            master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 4)
            typical_oracle_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)
            dcglp_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)
                            
            # solve mip for reference
            mip = Mip(data)
            assign_attributes!(mip.model, mip_solver_param)
            update_model!(mip, data)
            optimize!(mip.model)
            @assert termination_status(mip.model) == OPTIMAL
            mip_opt_val = objective_value(mip.model)

            x_opt = value.(mip.model[:x])
            t_opt = value.(mip.model[:t])
            @debug x_opt
            @debug t_opt

            @testset "Classical oracle" begin
                @testset "Seq" begin        
                    for strengthened in [true; false], add_benders_cuts_to_master in [true; false; 2], reuse_dcglp in [true; false], p in [1.0; Inf], disjunctive_cut_append_rule in [NoDisjunctiveCuts(); AllDisjunctiveCuts(); DisjunctiveCutsSmallerIndices()]
                    # for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [Inf], disjunctive_cut_append_rule in [NoDisjunctiveCuts()]
                        @info "solving UFLP p$i - disjunctive oracle/classical - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule"
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule" begin
                            master = Master(data; solver_param = master_solver_param)
                            update_model!(master, data)

                            typical_oracles = [ClassicalOracle(data; solver_param = typical_oracle_solver_param); ClassicalOracle(data; solver_param = typical_oracle_solver_param)] # for kappa & nu
                            for k=1:2
                                update_model!(typical_oracles[k], data)
                            end

                            disjunctive_oracle = DisjunctiveOracle(data, typical_oracles; 
                                                                   solver_param = dcglp_solver_param,
                                                                   param = dcglp_param) 
                            oracle_param = DisjunctiveOracleParam(norm = LpNorm(p), 
                                                                    split_index_selection_rule = RandomFractional(),
                                                                    disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
                                                                    strengthened=strengthened, 
                                                                    add_benders_cuts_to_master=add_benders_cuts_to_master, 
                                                                    fraction_of_benders_cuts_to_master = 1.0, 
                                                                    reuse_dcglp=reuse_dcglp)
                            set_parameter!(disjunctive_oracle, oracle_param)
                            update_model!(disjunctive_oracle, data)

                            env = BendersSeq(data, master, disjunctive_oracle; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            # if env.log.termination_status == Optimal()
                            if !isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                                t_opt_ = [sum(t_opt)]
                                @error "Failed ***** mip_opt_val = $(mip_opt_val) vs BD_obj_val = $(env.obj_value)" 
                                # optimize!(master.model)
                                opt_sol = Dict{VariableRef, Float64}()
                                for i = 1:data.dim_x
                                    opt_sol[master.model[:x][i]] = x_opt[i]
                                end
                                for i = 1:data.dim_t
                                    opt_sol[master.model[:t][i]] = t_opt_[i]
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
                # SeqInOut for DisjunctiveOracle is actually unnecessary extra, but it works though
                @testset "SeqInOut" begin
                    # for strengthened in [false], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [Inf], disjunctive_cut_append_rule in [NoDisjunctiveCuts()]
                    for strengthened in [true; false], add_benders_cuts_to_master in [true; false; 2], reuse_dcglp in [true; false], p in [1.0; Inf], disjunctive_cut_append_rule in [NoDisjunctiveCuts(); AllDisjunctiveCuts(); DisjunctiveCutsSmallerIndices()]
                        @info "solving UFLP p$i - disjunctive oracle/classical - SeqInOut... - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule"
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule" begin
                            master = Master(data; solver_param = master_solver_param)
                            update_model!(master, data)

                            typical_oracles = [ClassicalOracle(data; solver_param = typical_oracle_solver_param); ClassicalOracle(data; solver_param = typical_oracle_solver_param)] # for kappa & nu
                            for k=1:2
                                update_model!(typical_oracles[k], data)
                            end

                            disjunctive_oracle = DisjunctiveOracle(data, typical_oracles; 
                                                                   solver_param = dcglp_solver_param,
                                                                   param = dcglp_param) 
                            oracle_param = DisjunctiveOracleParam(norm = LpNorm(p), 
                                                                    split_index_selection_rule = RandomFractional(),
                                                                    disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
                                                                    strengthened=strengthened, 
                                                                    add_benders_cuts_to_master=add_benders_cuts_to_master, 
                                                                    fraction_of_benders_cuts_to_master = 1.0, 
                                                                    reuse_dcglp=reuse_dcglp)
                            set_parameter!(disjunctive_oracle, oracle_param)
                            update_model!(disjunctive_oracle, data)

                            env = BendersSeqInOut(data, master, disjunctive_oracle; param = benders_inout_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            if !isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                                t_opt_ = [sum(t_opt)]
                                @error "Failed ***** mip_opt_val = $(mip_opt_val) vs BD_obj_val = $(env.obj_value)" 
                                # optimize!(master.model)
                                opt_sol = Dict{VariableRef, Float64}()
                                for i = 1:data.dim_x
                                    opt_sol[master.model[:x][i]] = x_opt[i]
                                end
                                for i = 1:data.dim_t
                                    opt_sol[master.model[:t][i]] = t_opt_[i]
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

            # initialize dim_x, dim_t, c_x, c_t
            dim_x = problem.n_facilities
            c_x = problem.fixed_costs
            dim_t = problem.n_customers # knapsack cut
            c_t = ones(dim_t)
            
            data = Data(dim_x, dim_t, problem, c_x, c_t)
            @assert dim_x == length(data.c_x)
            @assert dim_t == length(data.c_t)

            @testset "fat knapsack oracle" begin
                # fat-knapsack-based disjunctive cut has sparse gamma_t, so adding only disjunctive cut does not improve lower bound, setting add_benders_cuts_to_master = true
                @testset "Seq" begin
                    for strengthened in [true; false], add_benders_cuts_to_master in [true; 2], reuse_dcglp in [true; false], p in [1.0; Inf], disjunctive_cut_append_rule in [NoDisjunctiveCuts(); AllDisjunctiveCuts(); DisjunctiveCutsSmallerIndices()]
                    # for strengthened in [false], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [Inf], disjunctive_cut_append_rule in [NoDisjunctiveCuts()]
                        @info "solving UFLP p$i - disjunctive oracle/fatKnapsack -  Seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule"
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule" begin
                            master = Master(data; solver_param = master_solver_param)
                            update_model!(master, data)

                            # model-free knapsack-based cuts
                            typical_oracles = [UFLKnapsackOracle(data); UFLKnapsackOracle(data)] # for kappa & nu

                            disjunctive_oracle = DisjunctiveOracle(data, typical_oracles; 
                                                                   solver_param = dcglp_solver_param,
                                                                   param = dcglp_param) 
                            oracle_param = DisjunctiveOracleParam(norm = LpNorm(p), 
                                                                    split_index_selection_rule = RandomFractional(),
                                                                    disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
                                                                    strengthened=strengthened, 
                                                                    add_benders_cuts_to_master=add_benders_cuts_to_master, 
                                                                    fraction_of_benders_cuts_to_master = 1.0, 
                                                                    reuse_dcglp=reuse_dcglp)
                            # norm is used in the initialization.
                            set_parameter!(disjunctive_oracle, oracle_param)
                            update_model!(disjunctive_oracle, data)
                            
                            env = BendersSeq(data, master, disjunctive_oracle; param = benders_param)
                            log = solve!(env)
                            @test env.termination_status == Optimal()
                            # if env.log.termination_status == Optimal()
                            if !isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                                @error "Failed ***** mip_opt_val = $(mip_opt_val) vs BD_obj_val = $(env.obj_value)" 
                                # optimize!(master.model)
                                opt_sol = Dict{VariableRef, Float64}()
                                for i = 1:data.dim_x
                                    opt_sol[master.model[:x][i]] = x_opt[i]
                                end
                                for i = 1:data.dim_t
                                    opt_sol[master.model[:t][i]] = t_opt[i]
                                end
                                
                                @info primal_feasibility_report(env.master.model, opt_sol)
                                @info data.c_x' * x_opt + data.c_t' * t_opt
                                
                                for v in keys(opt_sol)
                                    fix(v, opt_sol[v]; force=true)
                                end
                                optimize!(env.master.model)
                                @info objective_value(env.master.model)
                            end
                            @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                            # elseif env.log.termination_status == TimeLimit()
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