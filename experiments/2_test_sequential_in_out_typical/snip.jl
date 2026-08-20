using BendersX
using CSV
using DataFrames
using Test
using JuMP
include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))

@testset verbose = true "SNIP Sequential Benders Tests" begin
    reference_path = normpath(joinpath(@__DIR__, "..", "reference_objectives", "snip.csv"))
    reference_df = DataFrame(CSV.File(reference_path))
    @assert nrow(reference_df) == length(unique(reference_df.instance_name)) "Duplicate SNIP reference objectives found in $(reference_path)"
    reference_objectives = Dict(String(row.instance_name) => Float64(row.objective_value) for row in eachrow(reference_df))
    for instance in [0], snipno in [0], budget in [30.0]
        @testset "instance $instance; snipno $snipno budget $budget" begin
            # Load problem data
            data = read_snip_data(instance, snipno, budget)

            # Loop parameters
            benders_inout_param = BendersSeqInOutParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = false,
                            stabilizing_x = ones(length(data.D)),
                            α = 0.9,
                            λ = 0.1
                        )

            instance_name = "instance=$(instance);snipno=$(snipno);budget=$(budget)"
            @assert haskey(reference_objectives, instance_name) "Missing SNIP reference objective for $(instance_name) in $(reference_path)"
            mip_opt_val = reference_objectives[instance_name]

            @testset "Classic oracle" begin     
                @info "solving SNIP instance-$instance snipno-$snipno budget-$budget - classical oracle - seq..."
                master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
                oracle = SeparableOracle(data, master, ClassicalOracle, data.num_scenarios; model = update_sub_model!, optimizer = optimizer)
                env = BendersSeqInOut(master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end 

            @testset "Pareto oracle" begin
                @info "solving SNIP instance-$instance snipno-$snipno budget-$budget - pareto oracle - seqInOut..."
                master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
                param = ParetoOracleParam(fill(1.0, length(data.D)))
                oracle = SeparableOracle(data, master, ParetoOracle, data.num_scenarios; model = update_sub_model!, sub_oracle_param = param, optimizer = optimizer)
                env = BendersSeqInOut(master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end

            @testset "Unified oracle" begin
                @info "solving SNIP instance-$instance snipno-$snipno budget-$budget - unified oracle - seqInOut..."
                master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
                oracle = SeparableOracle(data, master, UnifiedOracle, data.num_scenarios; model = update_sub_model!, sub_oracle_param = UnifiedOracleParam(), optimizer = optimizer)
                env = BendersSeqInOut(master, oracle; param = benders_inout_param)
                log = solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            end
        end
    end
end
