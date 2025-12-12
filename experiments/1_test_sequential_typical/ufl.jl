using BendersDecomposition
using Test
using JuMP
using CPLEX

@testset verbose = true "UFLP Sequential Benders Tests" begin
    instances = setdiff(1:71, [67])

    for i in instances
        @testset "Instance: p$i" begin
            
            # Load problem data
            data = read_uflp_benchmark_data("p$(i)")
            
            # Loop parameters
            benders_param = BendersSeqParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = false
                            )
            
            # Solve MIP for reference
            mip_model = Model()
            customize_mip_model!(mip_model, data)
            optimize!(mip_model)
            @assert termination_status(mip_model) == OPTIMAL
            mip_opt_val = objective_value(mip_model)

            @testset "Classic oracle" begin
                
                @info "solving UFLP p$i - classical oracle - seq..."
                
                master = Master(data; customize = customize_master_model!)
                oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end

            @testset "Knapsack oracle" begin
                
                function customize_master_model!(model::Model, data::UFLPData)
                    optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)
                    set_optimizer(model, optimizer)
                    @variable(model, x[1:data.n_facilities], Bin)
                    @variable(model, t[1:data.n_customers] >= -1e6)
                    @constraint(model, sum(x) >= 2)
                    @objective(model, Min, data.fixed_costs'* x + sum(t))
                    return (x = x, ), t
                end

                @testset "Fat version" begin
                    
                    @info "solving UFLP p$i - fat knapsack oracle - seq..."
                
                    master = Master(data; customize = customize_master_model!)
                    oracle = UFLKnapsackOracle(data) 
                    set_parameter!(oracle, "add_only_violated_cuts", true)

                    env = BendersSeq(master, oracle; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end

                # To test slim version, users can use # set_parameter!(oracle, "slim", true)
            end
        end
    end
end