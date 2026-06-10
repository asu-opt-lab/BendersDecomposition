using Test
using BendersX
using JuMP
using HiGHS
using MathOptInterface
using LinearAlgebra

const DNO_MOI = MathOptInterface

struct DisjunctiveNormTestData <: AbstractData
    costs::Matrix{Float64}
end

struct DirectionalVectorTTestData <: AbstractData end

struct DirectionalVectorTTestOracle <: BendersX.AbstractTypicalOracle end

function disjunctive_norm_optimizer()
    return optimizer_with_attributes(HiGHS.Optimizer, DNO_MOI.Silent() => true)
end

function disjunctive_norm_data()
    return DisjunctiveNormTestData([
        1.0 3.0
        3.0 1.0
    ])
end

function update_disjunctive_norm_master!(model::Model, data::DisjunctiveNormTestData)
    set_optimizer(model, disjunctive_norm_optimizer())
    n_facilities = size(data.costs, 1)

    @variable(model, x[1:n_facilities], Bin)
    @variable(model, t[1:1] >= -1.0e6)
    @objective(model, Min, 0.1 * sum(x) + t[1])
    @constraint(model, sum(x) >= 1.0)

    return (x = x,), t
end

function update_continuous_disjunctive_norm_master!(model::Model, data::DisjunctiveNormTestData)
    set_optimizer(model, disjunctive_norm_optimizer())
    n_facilities = size(data.costs, 1)

    @variable(model, 0 <= x[1:n_facilities] <= 1)
    @variable(model, t[1:1] >= -1.0e6)
    @objective(model, Min, 0.1 * sum(x) + t[1])
    @constraint(model, sum(x) >= 1.0)

    return (x = x,), t
end

function update_disjunctive_norm_sub!(model::Model, data::DisjunctiveNormTestData, scen_idx::Int; x)
    set_optimizer(model, disjunctive_norm_optimizer())
    n_facilities, n_customers = size(data.costs)

    @variable(model, y[1:n_facilities, 1:n_customers] >= 0)
    @objective(model, Min, sum(data.costs .* y))
    @constraint(model, demand[j in 1:n_customers], sum(y[:, j]) == 1.0)
    @constraint(model, facility_open[i in 1:n_facilities, j in 1:n_customers], y[i, j] <= x[i])

    return nothing
end

function build_disjunctive_norm_master(; continuous::Bool = false)
    data = disjunctive_norm_data()
    model_hook = continuous ? update_continuous_disjunctive_norm_master! : update_disjunctive_norm_master!
    return data, Master(data; model = model_hook, optimizer = disjunctive_norm_optimizer())
end

function build_typical_pair(data, master)
    return [
        ClassicalOracle(data, master; model = update_disjunctive_norm_sub!, optimizer = disjunctive_norm_optimizer()),
        ClassicalOracle(data, master; model = update_disjunctive_norm_sub!, optimizer = disjunctive_norm_optimizer()),
    ]
end

function disjunctive_norm_dcglp_param(; verbose::Bool = false)
    return DcglpParam(
        disjunctive_norm_optimizer();
        time_limit = 20.0,
        gap_tolerance = 1.0e-5,
        halt_limit = 3,
        iter_limit = 20,
        verbose = verbose,
    )
end

function update_directional_vector_t_master!(model::Model, ::DirectionalVectorTTestData)
    set_optimizer(model, disjunctive_norm_optimizer())

    @variable(model, x[1:2], Bin)
    @variable(model, t[1:2] >= -1.0e6)
    @objective(model, Min, sum(t))

    return (x = x,), t
end

function BendersX.generate_cuts(
    ::DirectionalVectorTTestOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64};
    tol_normalize::Float64 = 1.0,
    time_limit::Float64 = 3600.0,
)
    hyperplanes = BendersX.Hyperplane[]
    for j in eachindex(t_value)
        h = BendersX.Hyperplane(length(x_value), length(t_value))
        h.a_x[j] = -1.0
        h.a_t[j] = -1.0
        h.a_0 = 1.0
        push!(hyperplanes, h)
    end

    f_x = 1.0 .- x_value[1:length(t_value)]
    is_in_L = all(f_x .<= t_value .+ 1.0e-9)
    return is_in_L, hyperplanes, f_x
end

function run_direct_cut_smoke(oracle)
    is_in_L, hyperplanes, f_x = BendersX.generate_cuts(
        oracle,
        [0.5, 0.5],
        [0.0];
        time_limit = 20.0,
        throw_typical_cuts_for_errors = true,
    )

    @test is_in_L isa Bool
    @test hyperplanes isa Vector{BendersX.Hyperplane}
    @test f_x isa Vector{Float64}
    @test !isempty(hyperplanes)
    @test all(h -> BendersX.evaluate_violation(h, [0.5, 0.5], [0.0]) isa Number, hyperplanes)
end

@testset "Disjunctive norm oracles" begin
    @testset "evaluate_violation returns numeric value" begin
        h = BendersX.Hyperplane([1.0, -2.0], [-1.0], 0.5)
        @test BendersX.evaluate_violation(h, [2.0, 1.0], [0.25]) == 0.25
    end

    @testset "parameter validation" begin
        dcglp_param = disjunctive_norm_dcglp_param()
        @test_throws ArgumentError DistanceNormOracleParam(dcglp_param; add_benders_cuts_to_master = 3)
        @test_throws ArgumentError SimplexNormOracleParam(dcglp_param; fraction_of_benders_cuts_to_master = 0.0)
        @test_throws ArgumentError VerticalReversePolarOracleParam(dcglp_param; fraction_of_benders_cuts_to_master = 1.1)
        @test_throws ArgumentError DirectionalPolarOracleParam(dcglp_param, Float64[], [0.0])
    end

    @testset "constructor validation" begin
        data, master = build_disjunctive_norm_master()
        typical_oracles = build_typical_pair(data, master)
        dcglp_param = disjunctive_norm_dcglp_param()

        @test_throws ArgumentError DistanceNormOracle(master, typical_oracles[1:1], DistanceNormOracleParam(dcglp_param))
        @test_throws DimensionMismatch DirectionalPolarOracle(
            master,
            typical_oracles,
            DirectionalPolarOracleParam(dcglp_param, [0.25], [0.0]),
        )

        data_cont, continuous_master = build_disjunctive_norm_master(; continuous = true)
        continuous_typical_oracles = build_typical_pair(data_cont, continuous_master)
        @test_throws ArgumentError SimplexNormOracle(
            continuous_master,
            continuous_typical_oracles,
            SimplexNormOracleParam(dcglp_param),
        )
    end

    @testset "generic split oracle constructor" begin
        data, master = build_disjunctive_norm_master()
        param = DistanceNormOracleParam(disjunctive_norm_dcglp_param(); reuse_dcglp = false)
        oracle = SplitOracle(master, build_typical_pair(data, master), param)

        @test param isa DistanceNormOracleParam
        @test oracle isa SplitOracle
        @test oracle isa DistanceNormOracle
        @test oracle isa BendersX.AbstractSplitOracle
    end

    @testset "direct generate_cuts smoke" begin
        for (oracle_type, param) in [
            (DistanceNormOracle, DistanceNormOracleParam(disjunctive_norm_dcglp_param(); norm = LpNorm(Inf), reuse_dcglp = false)),
            (SimplexNormOracle, SimplexNormOracleParam(disjunctive_norm_dcglp_param(); reuse_dcglp = false)),
            (VerticalReversePolarOracle, VerticalReversePolarOracleParam(disjunctive_norm_dcglp_param(); reuse_dcglp = false)),
            (DirectionalPolarOracle, DirectionalPolarOracleParam(disjunctive_norm_dcglp_param(), [0.25, 0.25], [0.0]; reuse_dcglp = false)),
        ]
            data, master = build_disjunctive_norm_master()
            oracle = oracle_type(master, build_typical_pair(data, master), param)
            run_direct_cut_smoke(oracle)
        end
    end

    @testset "cut history and include flag" begin
        data, master = build_disjunctive_norm_master()
        param = VerticalReversePolarOracleParam(
            disjunctive_norm_dcglp_param();
            split_index_selection_rule = LargestFractional(),
            disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(),
            reuse_dcglp = false,
        )
        oracle = VerticalReversePolarOracle(master, build_typical_pair(data, master), param)

        _, hyperplanes, _ = BendersX.generate_cuts(
            oracle,
            [0.5, 0.5],
            [0.0];
            time_limit = 20.0,
            include_disjunctive_cuts_to_hyperplanes = false,
        )

        @test length(oracle.splits) == 1
        @test !isempty(oracle.disjunctive_cuts)
        @test all(cut -> !(cut in hyperplanes), oracle.disjunctive_cuts)
    end

    @testset "directional core point update" begin
        data, master = build_disjunctive_norm_master()
        oracle = DirectionalPolarOracle(
            master,
            build_typical_pair(data, master),
            DirectionalPolarOracleParam(disjunctive_norm_dcglp_param(), [0.25, 0.25], [0.0]),
        )

        set_core_point!(oracle, [0.2, 0.3], [0.1])
        @test oracle.param.core_point_x == [0.2, 0.3]
        @test oracle.param.core_point_t == [0.1]
        @test_throws DimensionMismatch set_core_point!(oracle, [0.1], [0.0])
    end

    @testset "directional lift cut normalization" begin
        data = DirectionalVectorTTestData()
        master = Master(data; model = update_directional_vector_t_master!, optimizer = disjunctive_norm_optimizer())
        param = DirectionalPolarOracleParam(
            disjunctive_norm_dcglp_param(),
            [0.5, 0.5],
            [0.75, 0.75];
            split_index_selection_rule = MostFractional(),
            disjunctive_cut_append_rule = AllDisjunctiveCuts(),
            add_benders_cuts_to_master = 2,
            reuse_dcglp = true,
            strengthened = false,
            lift = true,
        )
        oracle = DirectionalPolarOracle(master, [DirectionalVectorTTestOracle(), DirectionalVectorTTestOracle()], param)

        x_value = [0.5, 0.5]
        t_value = [0.0, 0.0]
        is_in_L, _, _ = BendersX.generate_cuts(oracle, x_value, t_value; time_limit = 20.0)

        @test !is_in_L
        @test !isempty(oracle.disjunctive_cuts)
        cut = last(oracle.disjunctive_cuts)
        direction_x = x_value .- oracle.param.core_point_x
        direction_t = t_value .- oracle.param.core_point_t
        @test isapprox(dot(cut.a_x, direction_x) + dot(cut.a_t, direction_t), 1.0; atol = 1.0e-6)
        @test BendersX.evaluate_violation(cut, x_value, t_value) > 0.0
    end

    @testset "BendersSeq solve smoke" begin
        data, master = build_disjunctive_norm_master()
        oracle = VerticalReversePolarOracle(
            master,
            build_typical_pair(data, master),
            VerticalReversePolarOracleParam(
                disjunctive_norm_dcglp_param();
                split_index_selection_rule = LargestFractional(),
                reuse_dcglp = false,
                add_benders_cuts_to_master = 2,
            ),
        )
        env = BendersSeq(master, oracle; param = BendersSeqParam(time_limit = 30.0, gap_tolerance = 1.0e-5, verbose = false))
        solve!(env)
        @test env.termination_status == Optimal()
        @test isfinite(env.obj_value)
    end
end
