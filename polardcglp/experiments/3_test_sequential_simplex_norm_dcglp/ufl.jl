using BendersX
using CSV
using DataFrames
using Test
using JuMP
using CPLEX

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))

@testset verbose = true "UFLP Sequential Simplex Norm DCGLP Tests -- MIP master" begin
    root_experiments_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "experiments"))
    reference_path = normpath(joinpath(root_experiments_dir, "reference_objectives", "uflp.csv"))
    reference_df = DataFrame(CSV.File(reference_path))
    @assert nrow(reference_df) == length(unique(reference_df.instance_name)) "Duplicate UFLP reference objectives found in $(reference_path)"
    reference_objectives = Dict(String(row.instance_name) => Float64(row.objective_value) for row in eachrow(reference_df))
    instances = setdiff(1:71, [67])

    function customize_master_knapsack!(model::Model, data::UFLPData)
        optimizer = optimizer_with_attributes(
            CPLEX.Optimizer,
            "CPXPARAM_Threads" => 7,
            "CPX_PARAM_EPINT" => 1e-9,
            "CPX_PARAM_EPRHS" => 1e-9,
            "CPX_PARAM_EPGAP" => 1e-6,
            MOI.Silent() => true,
        )
        set_optimizer(model, optimizer)
        @variable(model, x[1:data.n_facilities], Bin)
        @variable(model, t[1:data.n_customers] >= -1e6)
        @constraint(model, sum(x) >= 2)
        @objective(model, Min, data.fixed_costs' * x + sum(t))
        return (x = x,), t
    end

    for i in instances
        @testset "Instance: p$i" begin
            data = read_uflp_benchmark_data("p$i")

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
            @assert haskey(reference_objectives, instance_name) "Missing UFLP reference objective for $(instance_name) in $(reference_path)"
            mip_opt_val = reference_objectives[instance_name]

            @testset "Classical oracle" begin
                @testset "Seq" begin
                    @info "solving SimplexNormDCGLP UFLP p$i - classical - seq"

                    oracle_param = SimplexNormDCGLPParam(
                        dcglp_param;
                        split_index_selection_rule = LargestFractional(),
                        disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                        add_benders_cuts_to_master = false,
                        reuse_dcglp = false,
                        zero_tol = 1e-9,
                    )

                    master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)
                    typical_oracles = [
                        ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                        ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                    ]
                    disjunctive_oracle = SimplexNormDCGLPOracle(master, typical_oracles, oracle_param)
                    env = BendersSeq(master, disjunctive_oracle; param = benders_param)

                    solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                end
            end

            @testset "Fat knapsack oracle" begin
                @testset "Seq" begin
                    @info "solving SimplexNormDCGLP UFLP p$i - fat knapsack - seq"

                    oracle_param = SimplexNormDCGLPParam(
                        dcglp_param;
                        split_index_selection_rule = LargestFractional(),
                        disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                        add_benders_cuts_to_master = false,
                        reuse_dcglp = false,
                        zero_tol = 1e-9,
                    )

                    master = Master(data; customize = customize_master_knapsack!, optimizer = mip_optimizer)
                    typical_oracles = [UFLKnapsackOracle(data), UFLKnapsackOracle(data)]
                    disjunctive_oracle = SimplexNormDCGLPOracle(master, typical_oracles, oracle_param)
                    env = BendersSeq(master, disjunctive_oracle; param = benders_param)

                    solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                end
            end
        end
    end
end
