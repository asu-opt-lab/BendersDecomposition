using BendersX
using CSV
using DataFrames
using Test
using JuMP
using CPLEX
include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))

@testset verbose = true "UFLP Sequential Benders Tests" begin
    reference_path = normpath(joinpath(@__DIR__, "..", "reference_objectives", "uflp.csv"))
    reference_df = DataFrame(CSV.File(reference_path))
    @assert nrow(reference_df) == length(unique(reference_df.instance_name)) "Duplicate UFLP reference objectives found in $(reference_path)"
    reference_objectives = Dict(String(row.instance_name) => Float64(row.objective_value) for row in eachrow(reference_df))
    instances = 1:71

    # GBC-enabled subproblem customization (y[i,j] <= x[i] via GBC)
    function update_sub_gbc_model!(model::Model, data::UFLPData, scen_idx::Int; x)
        optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, MOI.Silent() => true)
        set_optimizer(model, optimizer)

        I, J = data.n_facilities, data.n_customers
        @variable(model, y[1:I, 1:J] >= 0)

        cost_demands = data.costs .* data.demands'
        @objective(model, Min, sum(cost_demands .* y))

        @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)

        # Return GBC tuple: y[i,j] <= x[i] for each j
        gbc_lhs = vec(y)
        gbc_rhs = [x[i] for j in 1:J for i in 1:I]  # j outer, i inner to match vec(y)
        gbc_sense = fill(UpperBound, I*J)
        return gbc_lhs, gbc_rhs, gbc_sense
    end

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

            instance_name = "p$i"
            @assert haskey(reference_objectives, instance_name) "Missing UFLP reference objective for $(instance_name) in $(reference_path)"
            mip_opt_val = reference_objectives[instance_name]

            @testset "Classic oracle" begin

                @info "solving UFLP p$i - classical oracle - seq..."

                master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
                oracle = ClassicalOracle(data, master; model = update_sub_model!, optimizer = optimizer)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end

            @testset "Unified oracle" begin
                @info "solving UFLP p$i - unified oracle - seq..."
                master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
                oracle = UnifiedOracle(data, master; model = update_sub_model!, optimizer = optimizer)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end

            @testset "Classic oracle with GBC" begin

                @info "solving UFLP p$i - classical oracle with GBC - seq..."

                master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
                oracle = ClassicalOracle(data, master; model = update_sub_gbc_model!, optimizer = optimizer)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end

            @testset "Knapsack oracle" begin

                function update_master_model!(model::Model, data::UFLPData)
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

                    master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
                    oracle = UFLKnapsackOracle(data)
                    set_parameter!(oracle, "add_only_violated_cuts", true)

                    env = BendersSeq(master, oracle; param = benders_param)
                    log = solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
                end

                # To test slim version, users can use # set_parameter!(oracle, "slim", true)
            end

            @testset "Pareto oracle" begin
                @info "solving UFLP p$i - pareto oracle - seq..."
                master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
                param = ParetoOracleParam(fill(1.0, data.n_facilities))
                oracle = ParetoOracle(data, master, param; model = update_sub_model!, optimizer = optimizer)
                env = BendersSeq(master, oracle; param = benders_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end
        end
    end
end
