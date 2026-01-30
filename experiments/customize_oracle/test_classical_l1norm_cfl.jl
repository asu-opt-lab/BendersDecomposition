"""
Test L1 Norm Feasibility Cuts in ClassicalOracle using CFLP

This test verifies that the L1 normalization for feasibility cuts in ClassicalOracle
works correctly by comparing against MIP reference solutions.

The L1 norm is now the default behavior for ClassicalOracle when the subproblem
is infeasible. When infeasibility is detected:
1. Model is dynamically transformed to L1 form (min z0)
2. Feasibility cut is extracted from dual multipliers
3. Model is restored to original form
"""

using BendersX
using Test
using JuMP
using CPLEX

function customize_master_model!(model::Model, data::CFLPData)
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, 
        "CPXPARAM_Threads" => 7, 
        "CPX_PARAM_EPINT" => 1e-9, 
        "CPX_PARAM_EPRHS" => 1e-9, 
        "CPX_PARAM_EPGAP" => 1e-6, 
        MOI.Silent() => true
    )
    set_optimizer(model, optimizer)
    
    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t >= -1e6)

    @objective(model, Min, data.fixed_costs' * x + t)
    # @constraint(model, capacity, sum(data.capacities[i] * x[i] for i in 1:data.n_facilities) >= sum(data.demands))

    return (x = x, ), t
end

function customize_sub_model!(model::Model, data::CFLPData, scen_idx::Int; x) 
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, 
        "CPXPARAM_Threads" => 7, 
        "CPX_PARAM_EPRHS" => 1e-9, 
        "CPX_PARAM_EPOPT" => 1e-9, 
        "CPX_PARAM_NUMERICALEMPHASIS" => 1, 
        MOI.Silent() => true
    )
    set_optimizer(model, optimizer)

    I, J = data.n_facilities, data.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    
    cost_demands = data.costs .* data.demands'
    @objective(model, Min, sum(cost_demands .* y))
    
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, capacity[i in 1:I], sum(data.demands[:] .* y[i,:]) <= data.capacities[i] * x[i])
    @constraint(model, linking[i in 1:I, j in 1:J], y[i, j] <= x[i])
    
    return nothing
end

function customize_sub_model_gbc!(model::Model, data::CFLPData, scen_idx::Int; x) 
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, 
        "CPXPARAM_Threads" => 7, 
        "CPX_PARAM_EPRHS" => 1e-9, 
        "CPX_PARAM_EPOPT" => 1e-9, 
        "CPX_PARAM_NUMERICALEMPHASIS" => 1, 
        MOI.Silent() => true
    )
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

@testset verbose = true "ClassicalOracle L1 Norm CFLP Tests" begin
    # Test on a subset of instances
    # instances = [1, 2, 3]
    instances = setdiff(1:71, [67])
    
    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data
            data = read_cflp_benchmark_data("p$i")
            
            # Benders parameters
            benders_param = BendersSeqParam(;
                time_limit = 200.0,
                gap_tolerance = 1e-6,
                verbose = false
            )

            # Solve MIP for reference optimal value
            mip_model = Model()
            customize_mip_model!(mip_model, data)
            optimize!(mip_model)
            @assert termination_status(mip_model) == OPTIMAL
            mip_opt_val = objective_value(mip_model)
            @info "MIP optimal value for p$i: $mip_opt_val"

            # Test ClassicalOracle (now uses L1 norm for feasibility cuts by default)
            @testset "ClassicalOracle (L1 norm default)" begin
                @info "solving CFLP p$i - ClassicalOracle with L1 norm..."
                master = Master(data; customize = customize_master_model!)
                oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                @info "ClassicalOracle: $(env.obj_value), iterations: $(size(log, 1))"
            end

            # Test ClassicalOracle with GBC
            @testset "ClassicalOracle with GBC (L1 norm)" begin
                @info "solving CFLP p$i - ClassicalOracle with GBC and L1 norm..."
                master = Master(data; customize = customize_master_model!)
                oracle = ClassicalOracle(data, master; customize = customize_sub_model_gbc!)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                @info "ClassicalOracle with GBC: $(env.obj_value), iterations: $(size(log, 1))"
            end

            # Compare with UnifiedOracle (should produce similar results)
            @testset "UnifiedOracle (comparison)" begin
                @info "solving CFLP p$i - UnifiedOracle for comparison..."
                master = Master(data; customize = customize_master_model!)
                oracle = UnifiedOracle(data, master; customize = customize_sub_model!)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                @info "UnifiedOracle: $(env.obj_value), iterations: $(size(log, 1))"
            end
        end
    end
end
