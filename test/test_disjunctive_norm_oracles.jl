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

struct ExtensionContractNormalization <: BendersX.AbstractNormalization end
struct MissingContractNormalization <: BendersX.AbstractNormalization end

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
    return (
        ClassicalOracle(data, master; model = update_disjunctive_norm_sub!, optimizer = disjunctive_norm_optimizer()),
        ClassicalOracle(data, master; model = update_disjunctive_norm_sub!, optimizer = disjunctive_norm_optimizer()),
    )
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

function BendersX.add_normalization_constraint!(
    normalization::ExtensionContractNormalization,
    master::AbstractMaster,
    dcglp::Model,
    tau::VariableRef,
    sx::AbstractVector{VariableRef},
    st::AbstractVector{VariableRef},
)
    var_vec = [tau; sx; st]
    @constraint(dcglp, con_extension_norm, var_vec in DNO_MOI.NormInfinityCone(length(var_vec)))
end

function BendersX.update_dcglp_upper_bound_and_gap!(
    ::ExtensionContractNormalization,
    state::BendersX.DcglpState,
    log::BendersX.DcglpLog,
    t_value::Vector{Float64},
)
    BendersX.update_upper_bound_and_gap!(
        state,
        log,
        (t1, t2) -> LinearAlgebra.norm([state.values[:sx]; t1 .+ t2 .- t_value], Inf),
    )
end

function BendersX.disjunctive_cut_normalization_value(
    normalization::ExtensionContractNormalization,
    dcglp::Model,
    gamma_x::Vector{Float64},
    gamma_t::Vector{Float64},
)
    return max(1.0, LinearAlgebra.norm(vcat(gamma_x, gamma_t), 1.0))
end

function run_direct_cut_smoke(oracle)
    is_in_L, hyperplanes, f_x = BendersX.generate_cuts(
        oracle,
        [0.5, 0.5],
        [0.0];
        time_limit = 20.0
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
        @test_throws ArgumentError SplitOracleParam(LpDistanceNormalization(); dcglp_param = dcglp_param, add_benders_cuts_to_master = 3)
        @test_throws ArgumentError SplitOracleParam(ReversePolarNormalization(); dcglp_param = dcglp_param, fraction_of_benders_cuts_to_master = 1.1)
        @test_throws ArgumentError ReversePolarNormalization(Float64[], [0.0])
        @test_throws ArgumentError ReversePolarNormalization(; core_point_x = [0.0])
        @test_throws ArgumentError ReversePolarNormalization(; core_direction_x = [0.0])
        @test_throws ArgumentError ReversePolarNormalization(; core_direction_x = [0.0], core_direction_t = [0.0])
        @test_throws ArgumentError ReversePolarNormalization(;
            core_point_x = [0.0],
            core_point_t = [0.0],
            core_direction_x = [0.0],
            core_direction_t = [1.0],
        )
    end

    @testset "normalization object constructor" begin
        default_dcglp_param = DcglpParam()
        @test default_dcglp_param isa DcglpParam

        default_reverse_polar = ReversePolarNormalization()
        @test default_reverse_polar.core_point_x === nothing
        @test default_reverse_polar.core_point_t === nothing
        @test default_reverse_polar.core_direction_x === nothing
        @test default_reverse_polar.core_direction_t === nothing
        @test !default_reverse_polar.use_core_point

        disjunctive_norm_param = LpDistanceNormalization(1.0)
        param = SplitOracleParam(
            disjunctive_norm_param;
            dcglp_param = disjunctive_norm_dcglp_param(),
            reuse_dcglp = false,
        )
        @test param isa SplitOracleParam
        @test param.normalization isa LpDistanceNormalization
        @test param.normalization.norm_p == 1.0
        @test !param.reuse_dcglp
        @test LpDistanceNormalization(1).norm_p == 1.0
        @test LpDistanceNormalization(Inf).norm_p == Inf
        @test_throws ArgumentError LpDistanceNormalization(3.0)

        adjusted_param = SplitOracleParam(
            ReversePolarNormalization();
            dcglp_param = disjunctive_norm_dcglp_param()
        )
        @test adjusted_param.normalization isa ReversePolarNormalization

        directional_param = SplitOracleParam(
            ReversePolarNormalization([0.25, 0.25], [0.0]);
            dcglp_param = disjunctive_norm_dcglp_param(),
        )
        @test directional_param isa SplitOracleParam
        @test directional_param.normalization isa ReversePolarNormalization
        @test directional_param.normalization.core_point_x == [0.25, 0.25]
        @test directional_param.normalization.use_core_point

        fixed_direction_param = SplitOracleParam(
            ReversePolarNormalization(; core_direction_x = [0, 0], core_direction_t = [1]);
            dcglp_param = disjunctive_norm_dcglp_param(),
        )
        @test fixed_direction_param isa SplitOracleParam
        @test fixed_direction_param.normalization isa ReversePolarNormalization
        @test fixed_direction_param.normalization.core_point_x === nothing
        @test fixed_direction_param.normalization.core_direction_x == [0.0, 0.0]
        @test fixed_direction_param.normalization.core_direction_t == [1.0]
        @test !fixed_direction_param.normalization.use_core_point
    end

    @testset "constructor validation" begin
        data, master = build_disjunctive_norm_master()
        typical_oracles = build_typical_pair(data, master)
        dcglp_param = disjunctive_norm_dcglp_param()

        @test_throws MethodError SplitOracle(
            master,
            collect(typical_oracles);
            param = SplitOracleParam(LpDistanceNormalization(); dcglp_param = dcglp_param),
        )
    end

    @testset "reverse polar initialization validates dimensions" begin
        data, master = build_disjunctive_norm_master()
        typical_oracles = build_typical_pair(data, master)
        dcglp_param = disjunctive_norm_dcglp_param()

        default_normalization = ReversePolarNormalization()
        default_oracle = SplitOracle(
            master,
            typical_oracles;
            param = SplitOracleParam(default_normalization; dcglp_param = dcglp_param),
        )
        @test default_oracle.param.normalization.core_direction_x == zeros(master.dim_x)
        @test default_oracle.param.normalization.core_direction_t == ones(master.dim_t)

        @test_throws DimensionMismatch SplitOracle(
            master,
            typical_oracles;
            param = SplitOracleParam(
                ReversePolarNormalization([0.25], [0.0]);
                dcglp_param = dcglp_param,
            ),
        )
        @test_throws DimensionMismatch SplitOracle(
            master,
            typical_oracles;
            param = SplitOracleParam(
                ReversePolarNormalization([0.25, 0.25], [0.0, 0.0]);
                dcglp_param = dcglp_param,
            ),
        )
        @test_throws DimensionMismatch SplitOracle(
            master,
            typical_oracles;
            param = SplitOracleParam(
                ReversePolarNormalization(; core_direction_x = [0.0], core_direction_t = [1.0]);
                dcglp_param = dcglp_param,
            ),
        )
        @test_throws DimensionMismatch SplitOracle(
            master,
            typical_oracles;
            param = SplitOracleParam(
                ReversePolarNormalization(; core_direction_x = [0.0, 0.0], core_direction_t = [1.0, 1.0]);
                dcglp_param = dcglp_param,
            ),
        )
    end

    @testset "generic split oracle constructor" begin
        data, master = build_disjunctive_norm_master()
        default_oracle = SplitOracle(master, build_typical_pair(data, master))
        @test default_oracle isa SplitOracle
        @test default_oracle.param.normalization isa LpDistanceNormalization
        @test default_oracle.param.normalization.norm_p == Inf

        oracle = SplitOracle(
            master,
            build_typical_pair(data, master);
            param = SplitOracleParam(
                LpDistanceNormalization();
                dcglp_param = disjunctive_norm_dcglp_param(),
                reuse_dcglp = false,
            ),
        )

        @test oracle isa SplitOracle
        @test oracle.param.normalization isa LpDistanceNormalization
        @test oracle isa BendersX.AbstractDisjunctiveOracle
        @test !oracle.param.reuse_dcglp
    end

    @testset "direct generate_cuts smoke" begin
        for normalization in [
            LpDistanceNormalization(Inf),
            ReversePolarNormalization(),
            ReversePolarNormalization([0.25, 0.25], [0.0]),
            ReversePolarNormalization(; core_direction_x = [0.0, 0.0], core_direction_t = [1.0]),
        ]
            data, master = build_disjunctive_norm_master()
            oracle = SplitOracle(
                master,
                build_typical_pair(data, master);
                param = SplitOracleParam(
                    normalization;
                    dcglp_param = disjunctive_norm_dcglp_param(),
                    reuse_dcglp = false,
                ),
            )
            run_direct_cut_smoke(oracle)
        end
    end

    @testset "normalization extension defaults" begin
        
        data, master = build_disjunctive_norm_master()

        model = Model()
        @variable(model, tau)
        @variable(model, sx[1:1])
        @variable(model, st[1:1])
        @test_throws BendersX.UnimplementedInterfaceException BendersX.add_normalization_constraint!(
            MissingContractNormalization(),
            master,
            model,
            tau,
            sx,
            st,
        )
        @test_throws BendersX.UnimplementedInterfaceException BendersX.update_dcglp_upper_bound_and_gap!(
            MissingContractNormalization(),
            BendersX.DcglpState(),
            BendersX.DcglpLog(),
            [0.0],
            [0.0],
        )
        @test_throws BendersX.UnimplementedInterfaceException BendersX.disjunctive_cut_normalization_value(
            MissingContractNormalization(),
            model,
            [0.0],
            [0.0],
        )

        oracle = SplitOracle(
            master,
            build_typical_pair(data, master);
            param = SplitOracleParam(
                ExtensionContractNormalization();
                dcglp_param = disjunctive_norm_dcglp_param(),
                reuse_dcglp = false,
            ),
        )

        run_direct_cut_smoke(oracle)
    end

    @testset "cut history and include flag" begin
        data, master = build_disjunctive_norm_master()
        oracle = SplitOracle(
            master,
            build_typical_pair(data, master);
            param = SplitOracleParam(
                ReversePolarNormalization();
                dcglp_param = disjunctive_norm_dcglp_param(),
                split_index_selection_rule = LargestFractional(),
                disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(),
                reuse_dcglp = false,
            ),
        )

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

    @testset "directional lift cut normalization" begin
        data = DirectionalVectorTTestData()
        master = Master(data; model = update_directional_vector_t_master!, optimizer = disjunctive_norm_optimizer())
        oracle = SplitOracle(
            master,
            (DirectionalVectorTTestOracle(), DirectionalVectorTTestOracle());
            param = SplitOracleParam(
                ReversePolarNormalization([0.5, 0.5], [0.75, 0.75]);
                dcglp_param = disjunctive_norm_dcglp_param(),
                split_index_selection_rule = MostFractional(),
                disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                add_benders_cuts_to_master = 2,
                reuse_dcglp = true,
                strengthened = false,
                lift = true,
            ),
        )

        x_value = [0.5, 0.5]
        t_value = [0.0, 0.0]
        is_in_L, _, _ = BendersX.generate_cuts(oracle, x_value, t_value; time_limit = 20.0)

        @test !is_in_L
        @test !isempty(oracle.disjunctive_cuts)
        cut = last(oracle.disjunctive_cuts)
        direction_x = x_value .- oracle.param.normalization.core_point_x
        direction_t = t_value .- oracle.param.normalization.core_point_t
        @test isapprox(dot(cut.a_x, direction_x) + dot(cut.a_t, direction_t), 1.0; atol = 1.0e-6)
        @test BendersX.evaluate_violation(cut, x_value, t_value) > 0.0
    end

    @testset "fixed direction lift cut normalization" begin
        data = DirectionalVectorTTestData()
        master = Master(data; model = update_directional_vector_t_master!, optimizer = disjunctive_norm_optimizer())
        direction_x = [0.0, 0.0]
        direction_t = [0.75, 0.75]
        oracle = SplitOracle(
            master,
            (DirectionalVectorTTestOracle(), DirectionalVectorTTestOracle());
            param = SplitOracleParam(
                ReversePolarNormalization(; core_direction_x = direction_x, core_direction_t = direction_t);
                dcglp_param = disjunctive_norm_dcglp_param(),
                split_index_selection_rule = MostFractional(),
                disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                add_benders_cuts_to_master = 2,
                reuse_dcglp = true,
                strengthened = false,
                lift = true,
            ),
        )

        x_value = [0.5, 0.5]
        t_value = [0.0, 0.0]
        is_in_L, _, _ = BendersX.generate_cuts(oracle, x_value, t_value; time_limit = 20.0)

        @test !is_in_L
        @test !isempty(oracle.disjunctive_cuts)
        cut = last(oracle.disjunctive_cuts)
        cut_direction_x = .-oracle.param.normalization.core_direction_x
        cut_direction_t = .-oracle.param.normalization.core_direction_t
        @test isapprox(dot(cut.a_x, cut_direction_x) + dot(cut.a_t, cut_direction_t), 1.0; atol = 1.0e-6)
        @test BendersX.evaluate_violation(cut, x_value, t_value) > 0.0
    end

    @testset "BendersSeq solve smoke" begin
        data, master = build_disjunctive_norm_master()
        oracle = SplitOracle(
            master,
            build_typical_pair(data, master);
            param = SplitOracleParam(
                ReversePolarNormalization();
                dcglp_param = disjunctive_norm_dcglp_param(),
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
