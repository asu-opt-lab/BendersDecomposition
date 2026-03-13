using BendersX
import BendersX: customize_sub_model!, UFLPData, UFLKnapsackOracle, read_uflp_benchmark_data
using CSV
using DataFrames
using Test
using JuMP
using CPLEX

@testset verbose = true "UFLP Sequential In/Out Benders Tests" begin
    reference_path = normpath(joinpath(@__DIR__, "..", "reference_objectives", "uflp.csv"))
    reference_df = DataFrame(CSV.File(reference_path))
    @assert nrow(reference_df) == length(unique(reference_df.instance_name)) "Duplicate UFLP reference objectives found in $(reference_path)"
    reference_objectives = Dict(String(row.instance_name) => Float64(row.objective_value) for row in eachrow(reference_df))
    instances = setdiff(1:71, [67])

    # GBC-enabled subproblem customization (y[i,j] <= x[i] via GBC)
    function customize_sub_model_gbc!(model::Model, data::UFLPData, scen_idx::Int; x) 
        optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, MOI.Silent() => true)
        set_optimizer(model, optimizer)

        I, J = data.n_facilities, data.n_customers
        @variable(model, y[1:I, 1:J] >= 0)

        cost_demands = data.costs .* data.demands'
        @objective(model, Min, sum(cost_demands .* y))

        @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
        
        # Return GBC tuple: y[i,j] <= x[i] for each j
        gbc_lhs = vec(y)
        gbc_rhs = [x[i] for j in 1:J for i in 1:I]
        gbc_sense = fill(UpperBound, I*J)
        return gbc_lhs, gbc_rhs, gbc_sense
    end

    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data
            data = read_uflp_benchmark_data("p$(i)")
            
            # Loop parameters
            benders_inout_param = BendersSeqInOutParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = false,
                            stabilizing_x = ones(data.n_facilities),
                            α = 0.9,
                            λ = 0.1
                        )
                        
            instance_name = "p$i"
            @assert haskey(reference_objectives, instance_name) "Missing UFLP reference objective for $(instance_name) in $(reference_path)"
            mip_opt_val = reference_objectives[instance_name]

            @testset "Classic oracle" begin
                @info "solving UFLP p$i - classical oracle - seqInOut..."
                master = Master(data; customize = customize_master_model!)
                oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
                env = BendersSeqInOut(master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end

            @testset "Unified oracle" begin
                @info "solving UFLP p$i - unified oracle - seqInOut..."
                master = Master(data; customize = customize_master_model!)
                oracle = UnifiedOracle(data, master; customize = customize_sub_model!)
                env = BendersSeqInOut(master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end

            @testset "Classic oracle with GBC" begin
                @info "solving UFLP p$i - classical oracle with GBC - seqInOut..."
                master = Master(data; customize = customize_master_model!)
                oracle = ClassicalOracle(data, master; customize = customize_sub_model_gbc!)
                env = BendersSeqInOut(master, oracle; param = benders_inout_param)
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

                    env = BendersSeqInOut(master, oracle; param = benders_inout_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end

                # To test slim version, users can use # set_parameter!(oracle, "slim", true)
            end

            @testset "Pareto oracle" begin
                @info "solving UFLP p$i - pareto oracle - seqInOut..."
                master = Master(data; customize = customize_master_model!)
                param = ParetoOracleParam(fill(1.0, data.n_facilities))
                oracle = ParetoOracle(data, master, param; customize = customize_sub_model!)
                env = BendersSeqInOut(master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end
        end
    end
end
