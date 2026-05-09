using Test
using JuMP
using HiGHS
using BendersX
using MathOptInterface
using LinearAlgebra

isdefined(Main, :SimplexNormDCGLP) || include(joinpath(@__DIR__, "..", "src", "SimplexNormDCGLP.jl"))
using .SimplexNormDCGLP

const RMOI = MathOptInterface

struct ReversePolarData <: AbstractData
    open_cost::Vector{Float64}
    ship_cost::Matrix{Float64}
end

reversepolar_data() = ReversePolarData(
    [0.5, 0.6, 0.8],
    [1.0 3.0 2.0 4.0;
     2.0 1.0 2.5 2.0;
     3.0 2.0 1.0 1.5],
)

reversepolar_highs_optimizer() = optimizer_with_attributes(HiGHS.Optimizer, RMOI.Silent() => true)

function customize_master_reversepolar!(model::Model, data::ReversePolarData)
    set_optimizer(model, reversepolar_highs_optimizer())

    n_facilities = length(data.open_cost)
    @variable(model, x[1:n_facilities], Bin)
    @variable(model, t >= -1e6)
    @objective(model, Min, dot(data.open_cost, x) + t)
    @constraint(model, sum(x) >= 1)
    return (x = x,), t
end

function customize_sub_reversepolar!(model::Model, data::ReversePolarData, scen_idx::Int; x)
    set_optimizer(model, reversepolar_highs_optimizer())

    n_facilities, n_customers = size(data.ship_cost)
    @variable(model, y[1:n_facilities, 1:n_customers] >= 0)
    @objective(model, Min, sum(data.ship_cost .* y))
    @constraint(model, demand[j in 1:n_customers], sum(y[:, j]) == 1)
    @constraint(model, facility_open[i in 1:n_facilities, j in 1:n_customers], y[i, j] <= x[i])
    return nothing
end

function solve_reversepolar_reference(data::ReversePolarData)
    model = Model()
    set_optimizer(model, reversepolar_highs_optimizer())
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

function make_reversepolar_typical_oracles(data::ReversePolarData, master::Master)
    return [
        ClassicalOracle(data, master; customize = customize_sub_reversepolar!),
        ClassicalOracle(data, master; customize = customize_sub_reversepolar!),
    ]
end

function make_reversepolar_oracle(
    data::ReversePolarData,
    master::Master;
    split_rule = MostFractional(),
    append_rule = AllDisjunctiveCuts(),
    reuse_dcglp::Bool = true,
    verbose::Bool = false,
)
    dcglp_param = DcglpParam(
        reversepolar_highs_optimizer();
        time_limit = 10.0,
        gap_tolerance = 1e-8,
        halt_limit = 3,
        iter_limit = 25,
        verbose = verbose,
    )
    param = ReversePolarDCGLPParam(
        dcglp_param;
        split_index_selection_rule = split_rule,
        disjunctive_cut_append_rule = append_rule,
        add_benders_cuts_to_master = 2,
        fraction_of_benders_cuts_to_master = 1.0,
        reuse_dcglp = reuse_dcglp,
        zero_tol = 1e-9,
    )
    return ReversePolarDCGLPOracle(master, make_reversepolar_typical_oracles(data, master), param)
end

struct ReversePolarVectorTData <: AbstractData end

struct ReversePolarVectorTOracle <: BendersX.AbstractTypicalOracle end

function customize_master_reversepolar_vector_t!(model::Model, data::ReversePolarVectorTData)
    set_optimizer(model, reversepolar_highs_optimizer())
    @variable(model, x[1:2], Bin)
    @variable(model, t[1:2] >= -1e6)
    @objective(model, Min, sum(t))
    return (x = x,), t
end

function BendersX.generate_cuts(
    ::ReversePolarVectorTOracle,
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

@testset "ReversePolarDCGLP Environment" begin
    @testset "Constructor Validation" begin
        data = reversepolar_data()
        master = Master(data; customize = customize_master_reversepolar!)
        oracle = make_reversepolar_oracle(data, master)
        @test oracle isa ReversePolarDCGLPOracle
        @test length(oracle.dcglp[:tau]) == 1
    end

    @testset "Direct Vector-T Cut Generation" begin
        vt_data = ReversePolarVectorTData()
        vt_master = Master(vt_data; customize = customize_master_reversepolar_vector_t!)
        dcglp_param = DcglpParam(
            reversepolar_highs_optimizer();
            time_limit = 10.0,
            gap_tolerance = 1e-8,
            halt_limit = 3,
            iter_limit = 25,
            verbose = false,
        )
        param = ReversePolarDCGLPParam(
            dcglp_param;
            split_index_selection_rule = MostFractional(),
            disjunctive_cut_append_rule = AllDisjunctiveCuts(),
            add_benders_cuts_to_master = 2,
            reuse_dcglp = true,
            strengthened = false,
            zero_tol = 1e-9,
        )
        oracle = ReversePolarDCGLPOracle(vt_master, [ReversePolarVectorTOracle(), ReversePolarVectorTOracle()], param)

        x_value = zeros(2)
        t_value = zeros(2)
        is_in_L, _, _ = BendersX.generate_cuts(oracle, x_value, t_value; time_limit = 10.0)
        @test !is_in_L
        @test !isempty(oracle.disjunctiveCuts)

        cut = last(oracle.disjunctiveCuts)
        @test length(cut.a_t) == 2
        @test SimplexNormDCGLP.hyperplane_violation(cut, x_value, t_value) > 0.0
    end

    @testset "Fallback to Typical Cuts" begin
        vt_data = ReversePolarVectorTData()
        vt_master = Master(vt_data; customize = customize_master_reversepolar_vector_t!)
        dcglp_param = DcglpParam(
            reversepolar_highs_optimizer();
            time_limit = 10.0,
            gap_tolerance = 1e-8,
            halt_limit = 3,
            iter_limit = 25,
            verbose = false,
        )
        param = ReversePolarDCGLPParam(
            dcglp_param;
            split_index_selection_rule = MostFractional(),
            disjunctive_cut_append_rule = AllDisjunctiveCuts(),
            add_benders_cuts_to_master = 2,
            reuse_dcglp = true,
            strengthened = false,
            zero_tol = 1e-9,
        )
        oracle = ReversePolarDCGLPOracle(vt_master, [ReversePolarVectorTOracle(), ReversePolarVectorTOracle()], param)

        x_value = ones(2)
        t_value = zeros(2)
        is_in_L, hyperplanes, f_x = BendersX.generate_cuts(oracle, x_value, t_value; time_limit = 10.0)
        @test is_in_L
        @test isempty(oracle.disjunctiveCuts)
        @test length(hyperplanes) == 2
        @test all(isapprox.(f_x, 0.0; atol = 1e-9))
    end

    @testset "Verbose Smoke" begin
        vt_data = ReversePolarVectorTData()
        vt_master = Master(vt_data; customize = customize_master_reversepolar_vector_t!)
        dcglp_param = DcglpParam(
            reversepolar_highs_optimizer();
            time_limit = 10.0,
            gap_tolerance = 1e-8,
            halt_limit = 3,
            iter_limit = 25,
            verbose = true,
        )
        param = ReversePolarDCGLPParam(
            dcglp_param;
            split_index_selection_rule = MostFractional(),
            disjunctive_cut_append_rule = AllDisjunctiveCuts(),
            add_benders_cuts_to_master = 2,
            reuse_dcglp = true,
            strengthened = false,
            zero_tol = 1e-9,
        )
        oracle = ReversePolarDCGLPOracle(vt_master, [ReversePolarVectorTOracle(), ReversePolarVectorTOracle()], param)
        @test begin
            is_in_L, _, _ = BendersX.generate_cuts(oracle, zeros(2), zeros(2); time_limit = 10.0)
            !is_in_L
        end
    end

    @testset "BendersSeq Solve Path" begin
        data = reversepolar_data()
        reference_obj = solve_reversepolar_reference(data)

        master = Master(data; customize = customize_master_reversepolar!)
        oracle = make_reversepolar_oracle(data, master; reuse_dcglp = true)
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
