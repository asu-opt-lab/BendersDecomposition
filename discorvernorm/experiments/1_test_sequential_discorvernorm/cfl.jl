using BendersX
using CSV
using DataFrames
using Test
using JuMP
using CPLEX

include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))
include(normpath(joinpath(@__DIR__, "..", "discorvernorm_experiment_utils.jl")))
include(normpath(joinpath(@__DIR__, "..", "..", "src", "read_flcap.jl")))

@testset verbose = true "CFLP FLCAP Sequential Typical Tests" begin
    reference_path = normpath(joinpath(@__DIR__, "..", "reference_objectives", "cflp_flcap.csv"))
    reference_df = DataFrame(CSV.File(reference_path))
    @assert nrow(reference_df) == length(unique(reference_df.instance_name)) "Duplicate CFLP FLCAP reference objectives found in $(reference_path)"
    reference_objectives = Dict(String(row.instance_name) => Float64(row.objective_value) for row in eachrow(reference_df))
    instances = discorvernorm_flcap_instances()

    for instance_name in instances
        @testset "Instance: $instance_name" begin
            data = read_flcap_data(instance_name)

            benders_param = BendersSeqParam(
                time_limit = 800.0,
                gap_tolerance = 1e-6,
                verbose = false,
            )

            @assert haskey(reference_objectives, instance_name) "Missing CFLP FLCAP reference objective for $(instance_name) in $(reference_path)"
            mip_opt_val = reference_objectives[instance_name]

            @testset "Classical oracle" begin
                @info "solving typical CFLP FLCAP $instance_name - classical - seq"
                master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)
                oracle = ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer)
                env = BendersSeq(master, oracle; param = benders_param)

                solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol = 1.0, rtol = 1e-5)
            end

            @testset "Knapsack oracle" begin
                @info "solving typical CFLP FLCAP $instance_name - knapsack - seq"
                master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)
                oracle = CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer)
                env = BendersSeq(master, oracle; param = benders_param)

                solve!(env)
                @test env.termination_status == Optimal()
                @test isapprox(mip_opt_val, env.obj_value, atol = 1.0, rtol = 1e-5)
            end

        end
    end
end
