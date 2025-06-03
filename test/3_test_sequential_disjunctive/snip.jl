include("$(dirname(dirname(@__DIR__)))/example/snip/data_reader.jl")
include("$(dirname(dirname(@__DIR__)))/example/snip/model.jl")

@testset verbose = true "SNIP Sequential Benders Tests" begin
    for instance in [0], snipno in [0], budget in [30.0]
        @testset "instance $instance; snipno $snipno budget $budget" begin
            # Load problem data if necessary
            problem = read_snip_data(instance, snipno, budget)

            # initialize dim_x, dim_t, c_x, c_t
            dim_x = length(problem.D)
            dim_t = problem.num_scenarios
            c_x = zeros(dim_x)
            c_t = map(p -> p[3], problem.scenarios)
            data = Data(dim_x, dim_t, problem, c_x, c_t)
            @assert dim_x == length(data.c_x)
            @assert dim_t == length(data.c_t)

            # loop parameters
            benders_param = BendersSeqParam(;
                                            time_limit = 2000.0,
                                            gap_tolerance = 1e-6,
                                            verbose = true
                                            )
            dcglp_param = DcglpParam(;
                                    time_limit = 1000.0, 
                                    gap_tolerance = 1e-3, 
                                    halt_limit = 3, 
                                    iter_limit = 250,
                                    verbose = true
                                    )

            # solver parameters
            mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPXPARAM_Threads" => 4)
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

            @testset "Classic oracle" begin
                @testset "Seq" begin        
                    for strengthened in [true; false], add_benders_cuts_to_master in [true; 2], reuse_dcglp in [true; false], p in [1.0; Inf], disjunctive_cut_append_rule in [NoDisjunctiveCuts(); AllDisjunctiveCuts(); DisjunctiveCutsSmallerIndices()], adjust_t_to_fx in [true; false]  
                        @info "solving SNIP instance$instance snipno $snipno budget $budget - disjunctive oracle/classical - seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx"
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p dcut_append $disjunctive_cut_append_rule adjust_t_to_fx $adjust_t_to_fx" begin

                            master = Master(data; solver_param = master_solver_param)
                            update_model!(master, data)
                                     
                            typical_oracle_kappa = SeparableOracle(data, ClassicalOracle(), data.problem.num_scenarios; solver_param = typical_oracle_solver_param)
                            typical_oracle_nu = SeparableOracle(data, ClassicalOracle(), data.problem.num_scenarios; solver_param = typical_oracle_solver_param)
                            for j=1:typical_oracle_kappa.N
                                update_model!(typical_oracle_kappa.oracles[j], data, j)
                            end
                            for j=1:typical_oracle_nu.N
                                update_model!(typical_oracle_nu.oracles[j], data, j)
                            end
                            typical_oracles = [typical_oracle_kappa; typical_oracle_nu]

                            disjunctive_oracle = DisjunctiveOracle(data, typical_oracles; 
                                                                    solver_param = dcglp_solver_param,
                                                                    param = dcglp_param) 
                            oracle_param = DisjunctiveOracleParam(norm = LpNorm(p), 
                                                                    split_index_selection_rule = RandomFractional(),
                                                                    disjunctive_cut_append_rule = disjunctive_cut_append_rule, 
                                                                    strengthened=strengthened, 
                                                                    add_benders_cuts_to_master=add_benders_cuts_to_master, 
                                                                    fraction_of_benders_cuts_to_master = 1.0, 
                                                                    reuse_dcglp=reuse_dcglp,
                                                                    adjust_t_to_fx = adjust_t_to_fx)
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
            end 
            
        end
    end
end