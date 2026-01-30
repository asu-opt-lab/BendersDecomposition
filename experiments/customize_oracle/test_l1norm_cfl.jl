"""
Test L1NormOracle using CFLP (Capacitated Facility Location Problem)

This test compares L1NormOracle against ClassicalOracle to verify correctness.
"""

using BendersX
using Test
using JuMP
using CPLEX

# Include the L1NormOracle implementation
include("oracleL1Norm.jl")
function customize_master_model!(model::Model, data::CFLPData)
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)

    set_optimizer(model, optimizer)
    
    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t >= -1e6)

    @objective(model, Min, data.fixed_costs'* x + t)

    I = data.n_facilities
    # @constraint(model, capacity, sum(data.capacities[i] * x[i] for i in 1:I) >= sum(data.demands))

    return (x = x, ), t
end

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

@testset verbose = true "L1NormOracle CFLP Tests" begin
    # Test on a subset of instances for quick verification
    instances = [1]  # Can expand to more instances
    
    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data
            data = read_cflp_benchmark_data("p$i")
            
            # Benders parameters
            benders_param = BendersSeqParam(;
                time_limit = 200.0,
                gap_tolerance = 1e-6,
                verbose = true
            )

            # Solve MIP for reference optimal value
            mip_model = Model()
            customize_mip_model!(mip_model, data)
            optimize!(mip_model)
            @assert termination_status(mip_model) == OPTIMAL
            mip_opt_val = objective_value(mip_model)
            @info "MIP optimal value for p$i: $mip_opt_val"

            # Test Classical Oracle (baseline)
            # @testset "Classical Oracle (baseline)" begin
            #     @info "solving CFLP p$i - classical oracle..."
            #     master = Master(data; customize = customize_master_model!)
            #     oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
            #     env = BendersSeq(master, oracle; param = benders_param)
            #     log = solve!(env)
            #     @test env.termination_status == Optimal()
            #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            #     @info "Classical Oracle: $(env.obj_value), iterations: $(size(log, 1))"
            # end

            # Test L1Norm Oracle
            @testset "L1Norm Oracle" begin
                @info "solving CFLP p$i - L1Norm oracle..."
                master = Master(data; customize = customize_master_model!)
                oracle = L1NormOracle(data, master; customize = customize_sub_model_gbc!)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                @info "L1Norm Oracle: $(env.obj_value), iterations: $(size(log, 1))"
            end
        end
    end
end
