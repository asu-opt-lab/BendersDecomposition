using Test
using JuMP
using HiGHS
using BendersX
using MathOptInterface
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "DiscorverNormDCGLP.jl"))
import .DiscorverNormDCGLP

include(joinpath(@__DIR__, "..", "..", "polardcglp", "src", "PolarDCGLP.jl"))
import .PolarDCGLP

const DMOI = MathOptInterface

struct DiscorverSimplePolarData <: AbstractData
    open_cost::Vector{Float64}
    ship_cost::Matrix{Float64}
end

discorver_simple_polar_data() = DiscorverSimplePolarData(
    [0.5, 0.6, 0.8],
    [1.0 3.0 2.0 4.0;
     2.0 1.0 2.5 2.0;
     3.0 2.0 1.0 1.5],
)

discorver_highs_optimizer() = optimizer_with_attributes(HiGHS.Optimizer, DMOI.Silent() => true)

function customize_master_discorvernorm!(model::Model, data::DiscorverSimplePolarData)
    set_optimizer(model, discorver_highs_optimizer())

    n_facilities = length(data.open_cost)
    @variable(model, x[1:n_facilities], Bin)
    @variable(model, t >= -1e6)
    @objective(model, Min, dot(data.open_cost, x) + t)
    @constraint(model, sum(x) >= 1)
    return (x = x,), t
end

function customize_sub_discorvernorm!(model::Model, data::DiscorverSimplePolarData, scen_idx::Int; x)
    set_optimizer(model, discorver_highs_optimizer())

    n_facilities, n_customers = size(data.ship_cost)
    @variable(model, y[1:n_facilities, 1:n_customers] >= 0)
    @objective(model, Min, sum(data.ship_cost .* y))
    @constraint(model, demand[j in 1:n_customers], sum(y[:, j]) == 1)
    @constraint(model, facility_open[i in 1:n_facilities, j in 1:n_customers], y[i, j] <= x[i])
    return nothing
end

function solve_discorvernorm_reference(data::DiscorverSimplePolarData)
    model = Model()
    set_optimizer(model, discorver_highs_optimizer())
    n_facilities, n_customers = size(data.ship_cost)
    @variable(model, x[1:n_facilities], Bin)
    @variable(model, y[1:n_facilities, 1:n_customers] >= 0)
    @objective(model, Min, dot(data.open_cost, x) + sum(data.ship_cost .* y))
    @constraint(model, demand[j in 1:n_customers], sum(y[:, j]) == 1)
    @constraint(model, facility_open[i in 1:n_facilities, j in 1:n_customers], y[i, j] <= x[i])
    @constraint(model, sum(x) >= 1)
    optimize!(model)
    return objective_value(model)
end

function make_discorvernorm_typical_oracles(data::DiscorverSimplePolarData, master::Master)
    return [
        ClassicalOracle(data, master; customize = customize_sub_discorvernorm!),
        ClassicalOracle(data, master; customize = customize_sub_discorvernorm!),
    ]
end

function make_discorvernorm_oracle(
    data::DiscorverSimplePolarData,
    master::Master;
    support_points_x,
    support_points_t,
    split_rule = MostFractional(),
    append_rule = AllDisjunctiveCuts(),
    reuse_dcglp::Bool = true,
    strengthened::Bool = false,
    verbose::Bool = false,
)
    dcglp_param = DcglpParam(
        discorver_highs_optimizer();
        time_limit = 10.0,
        gap_tolerance = 1e-8,
        halt_limit = 3,
        iter_limit = 25,
        verbose = verbose,
    )
    param = DiscorverNormDCGLP.DiscorverNormDCGLPParam(
        dcglp_param,
        support_points_x,
        support_points_t;
        split_index_selection_rule = split_rule,
        disjunctive_cut_append_rule = append_rule,
        add_benders_cuts_to_master = 2,
        fraction_of_benders_cuts_to_master = 1.0,
        reuse_dcglp = reuse_dcglp,
        strengthened = strengthened,
        zero_tol = 1e-9,
    )
    return DiscorverNormDCGLP.DiscorverNormDCGLPOracle(master, make_discorvernorm_typical_oracles(data, master), param)
end

struct DiscorverVectorTData <: AbstractData end

struct DiscorverVectorTOracle <: BendersX.AbstractTypicalOracle end

function customize_master_discorvernorm_vector_t!(model::Model, data::DiscorverVectorTData)
    set_optimizer(model, discorver_highs_optimizer())
    @variable(model, x[1:2], Bin)
    @variable(model, t[1:2] >= -1e6)
    @objective(model, Min, sum(t))
    return (x = x,), t
end

function BendersX.generate_cuts(
    ::DiscorverVectorTOracle,
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
    is_in_L = all(f_x .<= t_value .+ 1e-9)
    return is_in_L, hyperplanes, f_x
end

function make_directional_reference_oracle(vt_master::Master, core_x::Vector{Float64}, core_t::Vector{Float64})
    dcglp_param = DcglpParam(
        discorver_highs_optimizer();
        time_limit = 10.0,
        gap_tolerance = 1e-8,
        halt_limit = 3,
        iter_limit = 25,
        verbose = false,
    )
    param = PolarDCGLP.DirectionalPolarDCGLPParam(
        dcglp_param,
        core_x,
        core_t;
        split_index_selection_rule = MostFractional(),
        disjunctive_cut_append_rule = AllDisjunctiveCuts(),
        add_benders_cuts_to_master = 2,
        reuse_dcglp = true,
        strengthened = false,
        zero_tol = 1e-9,
    )
    return PolarDCGLP.DirectionalPolarDCGLPOracle(vt_master, [DiscorverVectorTOracle(), DiscorverVectorTOracle()], param)
end

@testset "DiscorverNormDCGLP Environment" begin
    @testset "Constructor Validation" begin
        dcglp_param = DcglpParam(discorver_highs_optimizer(); verbose = false)
        @test_throws ArgumentError DiscorverNormDCGLP.DiscorverNormDCGLPParam(dcglp_param, zeros(2, 0), zeros(1, 0))
        @test_throws DimensionMismatch DiscorverNormDCGLP.DiscorverNormDCGLPParam(dcglp_param, zeros(2, 2), zeros(1, 1))

        data = discorver_simple_polar_data()
        master = Master(data; customize = customize_master_discorvernorm!)
        valid_support_x = hcat(fill(0.6, master.dim_x), ones(master.dim_x))
        valid_support_t = hcat([20.0], [20.0])
        oracle = make_discorvernorm_oracle(data, master; support_points_x = valid_support_x, support_points_t = valid_support_t)
        @test oracle isa DiscorverNormDCGLP.DiscorverNormDCGLPOracle
        @test haskey(oracle.dcglp, :tau_var)
        @test haskey(oracle.dcglp, :lambda_var)

        @test_throws DimensionMismatch DiscorverNormDCGLP.DiscorverNormDCGLPOracle(
            master,
            make_discorvernorm_typical_oracles(data, master),
            DiscorverNormDCGLP.DiscorverNormDCGLPParam(dcglp_param, zeros(master.dim_x - 1, 1), zeros(master.dim_t, 1)),
        )

        struct NonBinaryData <: AbstractData end

        function customize_master_nonbinary!(model::Model, data::NonBinaryData)
            set_optimizer(model, discorver_highs_optimizer())
            @variable(model, x[1:2] >= 0)
            @variable(model, t >= 0)
            @objective(model, Min, t)
            return (x = x,), t
        end

        nb_data = NonBinaryData()
        nb_master = Master(nb_data; customize = customize_master_nonbinary!)
        @test_throws ArgumentError DiscorverNormDCGLP.DiscorverNormDCGLPOracle(
            nb_master,
            [DiscorverVectorTOracle(), DiscorverVectorTOracle()],
            DiscorverNormDCGLP.DiscorverNormDCGLPParam(dcglp_param, zeros(2, 1), zeros(1, 1)),
        )
    end

    @testset "k=1 Regression to Directional" begin
        vt_data = DiscorverVectorTData()
        vt_master = Master(vt_data; customize = customize_master_discorvernorm_vector_t!)

        core_x = [0.5, 0.5]
        core_t = [0.75, 0.75]
        support_points_x = reshape(core_x, :, 1)
        support_points_t = reshape(core_t, :, 1)

        discorver_dcglp_param = DcglpParam(
            discorver_highs_optimizer();
            time_limit = 10.0,
            gap_tolerance = 1e-8,
            halt_limit = 3,
            iter_limit = 25,
            verbose = false,
        )
        discorver_param = DiscorverNormDCGLP.DiscorverNormDCGLPParam(
            discorver_dcglp_param,
            support_points_x,
            support_points_t;
            split_index_selection_rule = MostFractional(),
            disjunctive_cut_append_rule = AllDisjunctiveCuts(),
            add_benders_cuts_to_master = 2,
            reuse_dcglp = true,
            strengthened = false,
            zero_tol = 1e-9,
        )
        discorver_oracle = DiscorverNormDCGLP.DiscorverNormDCGLPOracle(vt_master, [DiscorverVectorTOracle(), DiscorverVectorTOracle()], discorver_param)
        directional_oracle = make_directional_reference_oracle(vt_master, core_x, core_t)

        x_value = zeros(2)
        t_value = zeros(2)
        is_in_L_discorver, _, _ = BendersX.generate_cuts(discorver_oracle, x_value, t_value; time_limit = 10.0)
        is_in_L_directional, _, _ = BendersX.generate_cuts(directional_oracle, x_value, t_value; time_limit = 10.0)
        @test !is_in_L_discorver
        @test !is_in_L_directional

        cut_discorver = last(discorver_oracle.disjunctiveCuts)
        cut_directional = last(directional_oracle.disjunctiveCuts)

        @test isapprox(cut_discorver.a_0, cut_directional.a_0; atol = 1e-6)
        @test isapprox(sort(abs.(Array(cut_discorver.a_x))), sort(abs.(Array(cut_directional.a_x))); atol = 1e-6)
        @test isapprox(sort(abs.(Array(cut_discorver.a_t))), sort(abs.(Array(cut_directional.a_t))); atol = 1e-6)
    end

    @testset "Multi-Anchor Direct Cut Generation" begin
        vt_data = DiscorverVectorTData()
        vt_master = Master(vt_data; customize = customize_master_discorvernorm_vector_t!)

        valid_anchor_x = [0.5, 0.5]
        valid_anchor_t = [0.75, 0.75]
        invalid_anchor_x = [0.0, 0.0]
        invalid_anchor_t = [0.75, 0.75]

        @test all((1.0 .- valid_anchor_x) .<= valid_anchor_t .+ 1e-9)
        @test !all((1.0 .- invalid_anchor_x) .<= invalid_anchor_t .+ 1e-9)

        dcglp_param = DcglpParam(
            discorver_highs_optimizer();
            time_limit = 10.0,
            gap_tolerance = 1e-8,
            halt_limit = 3,
            iter_limit = 25,
            verbose = false,
        )
        param = DiscorverNormDCGLP.DiscorverNormDCGLPParam(
            dcglp_param,
            hcat(invalid_anchor_x, valid_anchor_x),
            hcat(invalid_anchor_t, valid_anchor_t);
            split_index_selection_rule = MostFractional(),
            disjunctive_cut_append_rule = AllDisjunctiveCuts(),
            add_benders_cuts_to_master = 2,
            reuse_dcglp = true,
            strengthened = false,
            zero_tol = 1e-9,
        )
        oracle = DiscorverNormDCGLP.DiscorverNormDCGLPOracle(vt_master, [DiscorverVectorTOracle(), DiscorverVectorTOracle()], param)

        x_value = zeros(2)
        t_value = zeros(2)
        is_in_L, _, _ = BendersX.generate_cuts(oracle, x_value, t_value; time_limit = 10.0)
        @test !is_in_L
        @test !isempty(oracle.disjunctiveCuts)

        cut = last(oracle.disjunctiveCuts)
        support_diffs_x = hcat(x_value .- invalid_anchor_x, x_value .- valid_anchor_x)
        support_diffs_t = hcat(t_value .- invalid_anchor_t, t_value .- valid_anchor_t)
        @test isapprox(
            DiscorverNormDCGLP.compute_support_normalization_value(cut.a_x, cut.a_t, support_diffs_x, support_diffs_t),
            1.0;
            atol = 1e-6,
        )
        @test DiscorverNormDCGLP.hyperplane_violation(cut, x_value, t_value) > 0.0
    end

    @testset "BendersSeq Solve Path" begin
        data = discorver_simple_polar_data()
        reference_obj = solve_discorvernorm_reference(data)

        master = Master(data; customize = customize_master_discorvernorm!)
        support_points_x = hcat(fill(0.2, master.dim_x), ones(master.dim_x))
        support_points_t = hcat([20.0], [20.0])
        oracle = make_discorvernorm_oracle(
            data,
            master;
            support_points_x = support_points_x,
            support_points_t = support_points_t,
            reuse_dcglp = true,
            strengthened = false,
        )
        env = BendersSeq(
            master,
            oracle;
            param = BendersSeqParam(time_limit = 30.0, gap_tolerance = 1e-6, verbose = false),
        )
        solve!(env)
        @test env.termination_status == Optimal()
        @test isapprox(env.obj_value, reference_obj; atol = 1e-4)
    end
end
