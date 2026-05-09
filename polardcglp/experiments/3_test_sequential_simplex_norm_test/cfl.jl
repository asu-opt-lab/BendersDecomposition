using BendersX
using CSV
using DataFrames
using Test
using JuMP
using CPLEX

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))

@testset verbose = true "CFLP Sequential Simplex Norm Test Tests" begin
    root_experiments_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "experiments"))
    reference_path = normpath(joinpath(root_experiments_dir, "reference_objectives", "cflp.csv"))
    reference_df = DataFrame(CSV.File(reference_path))
    @assert nrow(reference_df) == length(unique(reference_df.instance_name)) "Duplicate CFLP reference objectives found in $(reference_path)"
    reference_objectives = Dict(String(row.instance_name) => Float64(row.objective_value) for row in eachrow(reference_df))
    instances = setdiff(1:71, [67])

    for i in instances
        @testset "Instance: p$i" begin
            data = read_cflp_benchmark_data("p$i")

            benders_param = BendersSeqParam(
                time_limit = 800.0,
                gap_tolerance = 1e-6,
                verbose = true,
            )
            dcglp_param = DcglpParam(
                dcglp_optimizer;
                time_limit = 1000.0,
                gap_tolerance = 1e-3,
                halt_limit = 3,
                iter_limit = 250,
                verbose = false,
            )

            instance_name = "p$i"
            @assert haskey(reference_objectives, instance_name) "Missing CFLP reference objective for $(instance_name) in $(reference_path)"
            mip_opt_val = reference_objectives[instance_name]

            @testset "Classical oracle" begin
                @testset "Seq" begin
                    @info "solving SimplexNormTestDCGLP CFLP p$i - classical - seq"

                    oracle_param = SimplexNormTestDCGLPParam(
                        dcglp_param;
                        split_index_selection_rule = RandomFractional(),
                        disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                        add_benders_cuts_to_master = false,
                        reuse_dcglp = false,
                        zero_tol = 1e-9,
                    )

                    master = Master(data; customize = customize_master_model!, optimizer = master_optimizer)
                    typical_oracles = [
                        ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = sub_optimizer),
                        ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = sub_optimizer),
                    ]
                    disjunctive_oracle = SimplexNormTestDCGLPOracle(master, typical_oracles, oracle_param)
                    env = BendersSeq(master, disjunctive_oracle; param = benders_param)

                    solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                end
            end

            @testset "Knapsack oracle" begin
                @testset "Seq" begin
                    @info "solving SimplexNormTestDCGLP CFLP p$i - knapsack - seq"

                    oracle_param = SimplexNormTestDCGLPParam(
                        dcglp_param;
                        split_index_selection_rule = RandomFractional(),
                        disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                        add_benders_cuts_to_master = false,
                        reuse_dcglp = false,
                        zero_tol = 1e-9,
                    )

                    master = Master(data; customize = customize_master_model!, optimizer = master_optimizer)
                    typical_oracles = [
                        CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = sub_optimizer),
                        CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = sub_optimizer),
                    ]
                    disjunctive_oracle = SimplexNormTestDCGLPOracle(master, typical_oracles, oracle_param)
                    env = BendersSeq(master, disjunctive_oracle; param = benders_param)

                    solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                end
            end
        end
    end
end
