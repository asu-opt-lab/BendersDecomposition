using BendersX
using CSV
using DataFrames
using Test
using JuMP
using CPLEX
using LinearAlgebra

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))

function build_cflp_core_point_x(data::CFLPData; optimizer = sub_optimizer)
    model = Model(optimizer)
    I = data.n_facilities

    @variable(model, delta >= 0.0)
    @variable(model, x[1:I])
    @objective(model, Max, delta)

    @constraint(model, [i in 1:I], x[i] >= delta)
    @constraint(model, [i in 1:I], x[i] <= 1.0 - delta)
    @constraint(model, sum(data.capacities[i] * x[i] for i in 1:I) >= sum(data.demands))
    @constraint(model, sum(x) >= 1.0)

    optimize!(model)

    if termination_status(model) != OPTIMAL
        throw(UnexpectedModelStatusException("Unable to build CFLP core point: termination status $(termination_status(model))."))
    end

    return value.(x), objective_value(model)
end

function build_cflp_ones_core_point_x(data::CFLPData)
    return ones(data.n_facilities), 0.0
end

function solve_cflp_recourse_value(data::CFLPData, x_value::Vector{Float64}; optimizer = sub_optimizer)
    model = Model(optimizer)
    I, J = data.n_facilities, data.n_customers
    @variable(model, y[1:I, 1:J] >= 0)

    cost_demands = data.costs .* data.demands'
    @objective(model, Min, sum(cost_demands .* y))
    @constraint(model, demand[j in 1:J], sum(y[:, j]) == 1)
    @constraint(model, facility_open[i in 1:I, j in 1:J], y[i, j] <= x_value[i])
    @constraint(model, capacity[i in 1:I], sum(data.demands .* y[i, :]) <= data.capacities[i] * x_value[i])

    optimize!(model)
    if termination_status(model) != OPTIMAL
        return Inf
    end
    return objective_value(model)
end

function build_directional_core_point(
    data::CFLPData;
    x_mode::String,
    t_margin_rel::Float64,
    t_margin_abs::Float64,
    optimizer = sub_optimizer,
)
    normalized_x_mode = lowercase(strip(x_mode))
    if normalized_x_mode in ("ones", "all_ones", "one")
        x_core, delta = build_cflp_ones_core_point_x(data)
    elseif normalized_x_mode in ("centered", "interior", "max_delta")
        x_core, delta = build_cflp_core_point_x(data; optimizer = optimizer)
    else
        throw(ArgumentError("Unsupported core_x_mode `$(x_mode)`. Use one of: ones, centered."))
    end

    recourse_value = solve_cflp_recourse_value(data, x_core; optimizer = optimizer)

    if !isfinite(recourse_value)
        @warn "Auto-generated CFLP core point was not subproblem-feasible; falling back to x = 1."
        x_core .= 1.0
        recourse_value = solve_cflp_recourse_value(data, x_core; optimizer = optimizer)
    end

    isfinite(recourse_value) || throw(UnexpectedModelStatusException("Unable to construct a finite core-point recourse value for CFLP."))

    t_margin = max(t_margin_abs, t_margin_rel * max(1.0, recourse_value))
    return x_core, [recourse_value + t_margin], delta, recourse_value
end

@testset verbose = true "CFLP Sequential Directional Polar Tests" begin
    root_experiments_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "experiments"))
    reference_path = normpath(joinpath(root_experiments_dir, "reference_objectives", "cflp.csv"))
    reference_df = DataFrame(CSV.File(reference_path))
    @assert nrow(reference_df) == length(unique(reference_df.instance_name)) "Duplicate CFLP reference objectives found in $(reference_path)"
    reference_objectives = Dict(String(row.instance_name) => Float64(row.objective_value) for row in eachrow(reference_df))
    instances = setdiff(1:71, [67])

    for i in instances
        @testset "Instance: p$i" begin
            data = read_cflp_benchmark_data("p$i")

            core_x, core_t, _, _ = build_directional_core_point(
                data;
                x_mode = "ones",
                t_margin_rel = 0.05,
                t_margin_abs = 1.0,
                optimizer = sub_optimizer,
            )

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
                    @info "solving DirectionalPolarDCGLP CFLP p$i - classical - seq"

                    oracle_param = DirectionalPolarDCGLPParam(
                        dcglp_param,
                        core_x,
                        core_t;
                        split_index_selection_rule = RandomFractional(),
                        disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                        add_benders_cuts_to_master = false,
                        reuse_dcglp = false,
                        zero_tol = 1e-9,
                    )

                    master = Master(data; customize = customize_master_model!, optimizer = master_optimizer)
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
            end

            @testset "Knapsack oracle" begin
                @testset "Seq" begin
                    @info "solving DirectionalPolarDCGLP CFLP p$i - knapsack - seq"

                    oracle_param = DirectionalPolarDCGLPParam(
                        dcglp_param,
                        core_x,
                        core_t;
                        split_index_selection_rule = RandomFractional(),
                        disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                        add_benders_cuts_to_master = false,
                        reuse_dcglp = false,
                        zero_tol = 1e-9,
                    )

                    master = Master(data; customize = customize_master_model!, optimizer = master_optimizer)
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
            end
        end
    end
end
