using BendersX
using Test
using JuMP
using CPLEX

@testset verbose = true "CFLP Sequential Benders Tests" begin
    instances = setdiff(1:71, [67])

    # GBC-enabled subproblem customization (y[i,j] <= x[i] via GBC)
    function customize_sub_model_gbc!(model::Model, data::CFLPData, scen_idx::Int; x) 
        optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, MOI.Silent() => true)
        set_optimizer(model, optimizer)

        I, J = data.n_facilities, data.n_customers
        @variable(model, y[1:I, 1:J] >= 0)
        
        cost_demands = data.costs .* data.demands'
        @objective(model, Min, sum(cost_demands .* y))
        
        @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
        @constraint(model, capacity[i in 1:I], sum(data.demands[:] .* y[i,:]) <= data.capacities[i] * x[i])
        
        # Return GBC tuple: y[i,j] <= x[i]
        gbc_lhs = vec(y)
        gbc_rhs = [x[i] for j in 1:J for i in 1:I]  # j outer, i inner to match vec(y)
        gbc_sense = fill(UpperBound, I*J)
        return gbc_lhs, gbc_rhs, gbc_sense
    end

    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data
            data = read_cflp_benchmark_data("p$i")
            
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
                @info "solving CFLP p$i - classical oracle - seq..."
                master = Master(data; customize = customize_master_model!)
                oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end 
            
            @testset "Knapsack oracle" begin
                @info "solving CFLP p$i - knapsack oracle - seq..."
                master = Master(data; customize = customize_master_model!)
                oracle = CFLKnapsackOracle(data, master; customize = customize_sub_model!)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end 

            @testset "Classic oracle with GBC" begin
                @info "solving CFLP p$i - classical oracle with GBC - seq..."
                master = Master(data; customize = customize_master_model!)
                oracle = ClassicalOracle(data, master; customize = customize_sub_model_gbc!)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end 

            @testset "Knapsack oracle with GBC" begin
                @info "solving CFLP p$i - knapsack oracle with GBC - seq..."
                master = Master(data; customize = customize_master_model!)
                oracle = CFLKnapsackOracle(data, master; customize = customize_sub_model_gbc!)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end

            @testset "Unified oracle" begin
                @info "solving CFLP p$i - unified oracle - seq..."
                master = Master(data; customize = customize_master_model!)
                oracle = UnifiedOracle(data, master; customize = customize_sub_model!)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end 
        end
    end
end