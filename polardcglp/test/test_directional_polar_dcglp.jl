using Test
using JuMP
using HiGHS
using BendersX
using MathOptInterface
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "DirectionalPolarDCGLP.jl"))
using .DirectionalPolarDCGLP

const DMOI = MathOptInterface

struct DirectionalSimplePolarData <: AbstractData
    open_cost::Vector{Float64}
    ship_cost::Matrix{Float64}
end

directional_simple_polar_data() = DirectionalSimplePolarData(
    [0.5, 0.6, 0.8],
    [1.0 3.0 2.0 4.0;
     2.0 1.0 2.5 2.0;
     3.0 2.0 1.0 1.5],
)

directional_highs_optimizer() = optimizer_with_attributes(HiGHS.Optimizer, DMOI.Silent() => true)

function customize_master_directional_polar!(model::Model, data::DirectionalSimplePolarData)
    set_optimizer(model, directional_highs_optimizer())

    n_facilities = length(data.open_cost)
    @variable(model, x[1:n_facilities], Bin)
    @variable(model, t >= -1e6)
    @objective(model, Min, dot(data.open_cost, x) + t)
    @constraint(model, sum(x) >= 1)
    return (x = x,), t
end

function customize_sub_directional_polar!(model::Model, data::DirectionalSimplePolarData, scen_idx::Int; x)
    set_optimizer(model, directional_highs_optimizer())

    n_facilities, n_customers = size(data.ship_cost)
    @variable(model, y[1:n_facilities, 1:n_customers] >= 0)
    @objective(model, Min, sum(data.ship_cost .* y))
    @constraint(model, demand[j in 1:n_customers], sum(y[:, j]) == 1)
    @constraint(model, facility_open[i in 1:n_facilities, j in 1:n_customers], y[i, j] <= x[i])
    return nothing
end

function solve_directional_reference(data::DirectionalSimplePolarData)
    model = Model()
    set_optimizer(model, directional_highs_optimizer())
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

function make_directional_typical_oracles(data::DirectionalSimplePolarData, master::Master)
    return [
        ClassicalOracle(data, master; customize = customize_sub_directional_polar!),
        ClassicalOracle(data, master; customize = customize_sub_directional_polar!),
    ]
end

function make_directional_polar_oracle(
    data::DirectionalSimplePolarData,
    master::Master;
    core_x = fill(0.6, master.dim_x),
    core_t = [10.0],
    split_rule = MostFractional(),
    append_rule = AllDisjunctiveCuts(),
    reuse_dcglp::Bool = true,
    strengthened::Bool = false,
    verbose::Bool = false,
)
    dcglp_param = DcglpParam(
        directional_highs_optimizer();
        time_limit = 10.0,
        gap_tolerance = 1e-8,
        halt_limit = 3,
        iter_limit = 25,
        verbose = verbose,
    )
    param = DirectionalPolarDCGLPParam(
        dcglp_param,
        collect(core_x),
        collect(core_t);
        split_index_selection_rule = split_rule,
        disjunctive_cut_append_rule = append_rule,
        add_benders_cuts_to_master = 2,
        fraction_of_benders_cuts_to_master = 1.0,
        reuse_dcglp = reuse_dcglp,
        strengthened = strengthened,
        zero_tol = 1e-9,
    )
    return DirectionalPolarDCGLPOracle(master, make_directional_typical_oracles(data, master), param)
end

struct DirectionalVectorTData <: AbstractData end

struct DirectionalVectorTOracle <: BendersX.AbstractTypicalOracle end

function customize_master_directional_vector_t!(model::Model, data::DirectionalVectorTData)
    set_optimizer(model, directional_highs_optimizer())
    @variable(model, x[1:2], Bin)
    @variable(model, t[1:2] >= -1e6)
    @objective(model, Min, sum(t))
    return (x = x,), t
end

function BendersX.generate_cuts(
    ::DirectionalVectorTOracle,
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

@testset "DirectionalPolarDCGLP Environment" begin
    @testset "Constructor Validation" begin
        data = directional_simple_polar_data()
        master = Master(data; customize = customize_master_directional_polar!)
        oracle = make_directional_polar_oracle(data, master)
        @test oracle isa DirectionalPolarDCGLPOracle
        @test haskey(oracle.dcglp, :tau_var)
        @test !haskey(oracle.dcglp, :lambda_var)

        dcglp_param = DcglpParam(directional_highs_optimizer(); verbose = false)
        @test_throws DimensionMismatch DirectionalPolarDCGLPOracle(
            master,
            make_directional_typical_oracles(data, master),
            DirectionalPolarDCGLPParam(dcglp_param, [0.5, 0.5], [10.0]),
        )
    end

    @testset "Direct Directional Cut Generation" begin
        vt_data = DirectionalVectorTData()
        vt_master = Master(vt_data; customize = customize_master_directional_vector_t!)
        dcglp_param = DcglpParam(
            directional_highs_optimizer();
            time_limit = 10.0,
            gap_tolerance = 1e-8,
            halt_limit = 3,
            iter_limit = 25,
            verbose = false,
        )
        param = DirectionalPolarDCGLPParam(
            dcglp_param,
            [0.5, 0.5],
            [0.75, 0.75];
            split_index_selection_rule = MostFractional(),
            disjunctive_cut_append_rule = AllDisjunctiveCuts(),
            add_benders_cuts_to_master = 2,
            reuse_dcglp = true,
            strengthened = false,
            zero_tol = 1e-9,
        )
        oracle = DirectionalPolarDCGLPOracle(vt_master, [DirectionalVectorTOracle(), DirectionalVectorTOracle()], param)

        x_value = zeros(2)
        t_value = zeros(2)
        is_in_L, _, _ = BendersX.generate_cuts(oracle, x_value, t_value; time_limit = 10.0)
        @test !is_in_L
        @test !isempty(oracle.disjunctiveCuts)

        cut = last(oracle.disjunctiveCuts)
        direction_x = x_value .- param.core_point_x
        direction_t = t_value .- param.core_point_t
        @test isapprox(dot(cut.a_x, direction_x) + dot(cut.a_t, direction_t), 1.0; atol = 1e-6)
        @test DirectionalPolarDCGLP.hyperplane_violation(cut, x_value, t_value) > 0.0
    end

    @testset "BendersSeq Solve Path" begin
        data = directional_simple_polar_data()
        reference_obj = solve_directional_reference(data)

        master = Master(data; customize = customize_master_directional_polar!)
        oracle = make_directional_polar_oracle(data, master; reuse_dcglp = true, strengthened = false)
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
