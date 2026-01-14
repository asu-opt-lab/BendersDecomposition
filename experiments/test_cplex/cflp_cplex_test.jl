# CFLP Test with CPLEX Solver
# This test verifies that the BendersX package works correctly with CPLEX

using BendersX
using Test
using JuMP
using MathOptInterface
const MOI = MathOptInterface

# Load CPLEX - this triggers the extension
using CPLEX

# Define CPLEX-specific optimizers
const MIP_OPTIMIZER = optimizer_with_attributes(
    CPLEX.Optimizer, 
    "CPXPARAM_Threads" => 4, 
    "CPX_PARAM_EPINT" => 1e-9, 
    "CPX_PARAM_EPRHS" => 1e-9, 
    "CPX_PARAM_EPGAP" => 1e-6, 
    MOI.Silent() => true
)

const LP_OPTIMIZER = optimizer_with_attributes(
    CPLEX.Optimizer, 
    "CPXPARAM_Threads" => 4, 
    "CPX_PARAM_EPRHS" => 1e-9, 
    "CPX_PARAM_EPOPT" => 1e-9, 
    "CPX_PARAM_NUMERICALEMPHASIS" => 1, 
    MOI.Silent() => true
)

# Custom model functions that use CPLEX
function customize_mip_model_cplex!(model::Model, data::CFLPData)
    set_optimizer(model, MIP_OPTIMIZER)
    
    I, J = data.n_facilities, data.n_customers
    @variable(model, x[1:I], Bin)
    @variable(model, y[1:I, 1:J] >= 0)
    @variable(model, t)
    
    cost_demands = data.costs .* data.demands'
    @objective(model, Min, data.fixed_costs'* x + t)

    @constraint(model, t >= sum(cost_demands .* y))
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.demands .* y[i,:]) <= data.capacities[i] * x[i])
    @constraint(model, capacity_total, sum(data.capacities[i] * x[i] for i in 1:I) >= sum(data.demands))
end

function customize_master_model_cplex!(model::Model, data::CFLPData)
    set_optimizer(model, MIP_OPTIMIZER)
    
    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t >= -1e6)

    @objective(model, Min, data.fixed_costs'* x + t)

    I = data.n_facilities
    @constraint(model, capacity, sum(data.capacities[i] * x[i] for i in 1:I) >= sum(data.demands))

    return (x = x, ), t
end

function customize_sub_model_cplex!(model::Model, data::CFLPData, scen_idx::Int; x)
    set_optimizer(model, LP_OPTIMIZER)

    I, J = data.n_facilities, data.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    cost_demands = data.costs .* data.demands'
    @objective(model, Min, sum(cost_demands .* y))
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.demands[:] .* y[i,:]) <= data.capacities[i] * x[i])
end

@testset verbose = true "CFLP with CPLEX" begin
    # Test with a single instance for quick verification
    instances = [1, 2, 3]

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
            customize_mip_model_cplex!(mip_model, data)
            optimize!(mip_model)
            @assert termination_status(mip_model) == OPTIMAL
            mip_opt_val = objective_value(mip_model)

            @testset "Classic oracle" begin
                @info "solving CFLP p$i with CPLEX - classical oracle..."
                master = Master(data; customize = customize_master_model_cplex!)
                oracle = ClassicalOracle(data, master; customize = customize_sub_model_cplex!)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end 
        end
    end
end


