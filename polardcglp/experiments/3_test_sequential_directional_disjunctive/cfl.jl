using BendersX
using CSV
using DataFrames
using Test
using JuMP
using CPLEX

isdefined(Main, :DirectionalPolarDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "DirectionalPolarDCGLP.jl")))
using .DirectionalPolarDCGLP

include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))

@testset verbose = true "CFLP Sequential Directional Disjunctive Tests" begin
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
                verbose = false,
            )
            dcglp_optimizer = optimizer_with_attributes(
                CPLEX.Optimizer,
                "CPXPARAM_Threads" => 7,
                "CPX_PARAM_EPRHS" => 1e-9,
                "CPX_PARAM_NUMERICALEMPHASIS" => 1,
                "CPX_PARAM_EPOPT" => 1e-9,
                MOI.Silent() => true,
            )
            dcglp_param = DcglpParam(
                dcglp_optimizer;
                time_limit = 1000.0,
                gap_tolerance = 1e-3,
                halt_limit = 3,
                iter_limit = 250,
                verbose = false,
            )
            benders_inout_param = BendersSeqInOutParam(
                time_limit = 800.0,
                gap_tolerance = 1e-6,
                verbose = false,
                stabilizing_x = ones(data.n_facilities),
                α = 0.9,
                λ = 0.1,
            )

            instance_name = "p$i"
            @assert haskey(reference_objectives, instance_name) "Missing CFLP reference objective for $(instance_name) in $(reference_path)"
            mip_opt_val = reference_objectives[instance_name]

            @testset "Classical oracle" begin
                @testset "Seq" begin
                    @info "solving DirectionalPolarDCGLP CFLP p$i - classical - seq"

                    oracle_param = DirectionalPolarDCGLPParam(
                        dcglp_param,
                        ones(data.n_facilities),
                        ones(1);
                        split_index_selection_rule = RandomFractional(),
                        disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                        add_benders_cuts_to_master = true,
                        fraction_of_benders_cuts_to_master = 1.0,
                        reuse_dcglp = true,

                        zero_tol = 1e-9,
                    )

                    master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)
                    typical_oracles = [
                        ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                        ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                    ]
                    disjunctive_oracle = DirectionalPolarDCGLPOracle(master, typical_oracles, oracle_param)
                    env = BendersSeq(master, disjunctive_oracle; param = benders_param)

                    solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                end

                @testset "SeqInOut" begin
                    @info "solving DirectionalPolarDCGLP CFLP p$i - classical - seqinout"

                    oracle_param = DirectionalPolarDCGLPParam(
                        dcglp_param,
                        ones(data.n_facilities),
                        ones(1);
                        split_index_selection_rule = RandomFractional(),
                        disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                        add_benders_cuts_to_master = true,
                        fraction_of_benders_cuts_to_master = 1.0,
                        reuse_dcglp = true,

                        zero_tol = 1e-9,
                    )

                    master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)
                    typical_oracles = [
                        ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                        ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                    ]
                    disjunctive_oracle = DirectionalPolarDCGLPOracle(master, typical_oracles, oracle_param)
                    env = BendersSeqInOut(master, disjunctive_oracle; param = benders_inout_param)

                    solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                end
            end

            @testset "Knapsack oracle" begin
                @testset "Seq" begin
                    @info "solving DirectionalPolarDCGLP CFLP p$i - knapsack - seq"

                    oracle_param = DirectionalPolarDCGLPParam(
                        dcglp_param,
                        ones(data.n_facilities),
                        ones(1);
                        split_index_selection_rule = RandomFractional(),
                        disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                        add_benders_cuts_to_master = true,
                        fraction_of_benders_cuts_to_master = 1.0,
                        reuse_dcglp = true,

                        zero_tol = 1e-9,
                    )

                    master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)
                    typical_oracles = [
                        CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                        CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                    ]
                    disjunctive_oracle = DirectionalPolarDCGLPOracle(master, typical_oracles, oracle_param)
                    env = BendersSeq(master, disjunctive_oracle; param = benders_param)

                    solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                end

                @testset "SeqInOut" begin
                    @info "solving DirectionalPolarDCGLP CFLP p$i - knapsack - seqinout"

                    oracle_param = DirectionalPolarDCGLPParam(
                        dcglp_param,
                        ones(data.n_facilities),
                        ones(1);
                        split_index_selection_rule = RandomFractional(),
                        disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                        add_benders_cuts_to_master = true,
                        fraction_of_benders_cuts_to_master = 1.0,
                        reuse_dcglp = true,

                        zero_tol = 1e-9,
                    )

                    master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)
                    typical_oracles = [
                        CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                        CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                    ]
                    disjunctive_oracle = DirectionalPolarDCGLPOracle(master, typical_oracles, oracle_param)
                    env = BendersSeqInOut(master, disjunctive_oracle; param = benders_inout_param)

                    solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                end
            end
        end
    end
end
