using Test
using JuMP
using HiGHS
using BendersX
using MathOptInterface
using LinearAlgebra
using Random
using SparseArrays

include(joinpath(@__DIR__, "..", "src", "PolarDCGLP.jl"))
using .PolarDCGLP

const MOI = MathOptInterface

struct SimplePolarData <: AbstractData
    open_cost::Vector{Float64}
    ship_cost::Matrix{Float64}
end

function simple_polar_data()
    return SimplePolarData(
        [0.5, 0.6, 0.8],
        [1.0 3.0 2.0 4.0;
         2.0 1.0 2.5 2.0;
         3.0 2.0 1.0 1.5],
    )
end

function highs_optimizer()
    return optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
end

function customize_master_polar!(model::Model, data::SimplePolarData)
    set_optimizer(model, highs_optimizer())

    I = length(data.open_cost)
    @variable(model, x[1:I], Bin)
    @variable(model, t >= -1e6)
    @objective(model, Min, dot(data.open_cost, x) + t)
    @constraint(model, sum(x) >= 1)
    return (x = x,), t
end

function customize_sub_polar!(model::Model, data::SimplePolarData, scen_idx::Int; x)
    set_optimizer(model, highs_optimizer())

    I, J = size(data.ship_cost)
    @variable(model, y[1:I, 1:J] >= 0)
    @objective(model, Min, sum(data.ship_cost .* y))
    @constraint(model, demand[j in 1:J], sum(y[:, j]) == 1)
    @constraint(model, facility_open[i in 1:I, j in 1:J], y[i, j] <= x[i])
    return nothing
end

function solve_reference(data::SimplePolarData)
    model = Model()
    set_optimizer(model, highs_optimizer())
    I, J = size(data.ship_cost)
    @variable(model, x[1:I], Bin)
    @variable(model, y[1:I, 1:J] >= 0)
    @objective(model, Min, dot(data.open_cost, x) + sum(data.ship_cost .* y))
    @constraint(model, demand[j in 1:J], sum(y[:, j]) == 1)
    @constraint(model, facility_open[i in 1:I, j in 1:J], y[i, j] <= x[i])
    @constraint(model, sum(x) >= 1)
    optimize!(model)
    return objective_value(model)
end

function make_typical_oracles(data::SimplePolarData, master::Master)
    return [
        ClassicalOracle(data, master; customize = customize_sub_polar!),
        ClassicalOracle(data, master; customize = customize_sub_polar!),
    ]
end

function make_polar_oracle(
    data::SimplePolarData,
    master::Master;
    split_rule = MostFractional(),
    append_rule = AllDisjunctiveCuts(),
    reuse_dcglp::Bool = true,
    adjust_t_to_fx::Bool = false,
    verbose::Bool = false,
)
    dcglp_param = DcglpParam(
        highs_optimizer();
        time_limit = 10.0,
        gap_tolerance = 1e-8,
        halt_limit = 3,
        iter_limit = 25,
        verbose = verbose,
    )
    param = PolarDCGLPParam(
        dcglp_param;
        split_index_selection_rule = split_rule,
        disjunctive_cut_append_rule = append_rule,
        add_benders_cuts_to_master = 2,
        fraction_of_benders_cuts_to_master = 1.0,
        reuse_dcglp = reuse_dcglp,
        adjust_t_to_fx = adjust_t_to_fx,
        zero_tol = 1e-9,
    )
    return PolarDCGLPOracle(master, make_typical_oracles(data, master), param)
end

struct VectorTOracle <: BendersX.AbstractTypicalOracle end

function BendersX.generate_cuts(
    ::VectorTOracle,
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
    is_in_L = sum(f_x) <= sum(t_value) + 1e-9
    return is_in_L, hyperplanes, f_x
end

@testset "PolarDCGLP Environment" begin
    data = simple_polar_data()
    reference_obj = solve_reference(data)

    @testset "Split Rules" begin
        x_value = [0.1, 0.49, 0.7, 1.0]

        phi, phi_0 = PolarDCGLP.split_phi_and_rhs(x_value, LargestFractional(); zero_tol = 1e-9)
        @test phi_0 == 0.0
        @test findfirst(x -> !iszero(x), phi) == 3

        phi, _ = PolarDCGLP.split_phi_and_rhs(x_value, MostFractional(); zero_tol = 1e-9)
        @test findfirst(x -> !iszero(x), phi) == 2

        Random.seed!(1)
        phi, _ = PolarDCGLP.split_phi_and_rhs(x_value, RandomFractional(); zero_tol = 1e-9)
        @test findfirst(x -> !iszero(x), phi) in (1, 2, 3)
    end

    @testset "Constructor Validation" begin
        master = Master(data; customize = customize_master_polar!)
        @test make_polar_oracle(data, master) isa PolarDCGLPOracle

        struct NonBinaryData <: AbstractData end

        function customize_master_nonbinary!(model::Model, data::NonBinaryData)
            set_optimizer(model, highs_optimizer())
            @variable(model, x[1:2] >= 0)
            @variable(model, t >= 0)
            @objective(model, Min, t)
            return (x = x,), t
        end

        nb_data = NonBinaryData()
        nb_master = Master(nb_data; customize = customize_master_nonbinary!)
        dcglp_param = DcglpParam(highs_optimizer(); verbose = false)
        param = PolarDCGLPParam(dcglp_param)
        typical = [ClassicalOracle(data, master; customize = customize_sub_polar!), ClassicalOracle(data, master; customize = customize_sub_polar!)]
        @test_throws ArgumentError PolarDCGLPOracle(nb_master, typical, param)

        struct VectorTData <: AbstractData end

        function customize_master_vector_t!(model::Model, data::VectorTData)
            set_optimizer(model, highs_optimizer())
            @variable(model, x[1:2], Bin)
            @variable(model, t[1:2] >= -1e6)
            @objective(model, Min, sum(t))
            return (x = x,), t
        end

        vt_data = VectorTData()
        vt_master = Master(vt_data; customize = customize_master_vector_t!)
        vt_typical = [VectorTOracle(), VectorTOracle()]
        @test PolarDCGLPOracle(vt_master, vt_typical, param) isa PolarDCGLPOracle
    end

    @testset "Disjunctive Cut Reuse Rules" begin
        master = Master(data; customize = customize_master_polar!)
        oracle = make_polar_oracle(data, master; append_rule = DisjunctiveCutsSmallerIndices())
        cut_1 = BendersX.Hyperplane([1.0, 0.0, 0.0], [-1.0], 0.2)
        cut_2 = BendersX.Hyperplane([0.0, 1.0, 0.0], [-1.0], 0.4)
        push!(oracle.disjunctiveCutsByIndex[1], cut_1)
        push!(oracle.disjunctiveCutsByIndex[2], cut_2)
        push!(oracle.disjunctiveCuts, cut_1, cut_2)
        push!(oracle.splits, (sparsevec([2], [1.0], master.dim_x), 0.0))
        PolarDCGLP.add_disjunctive_cuts!(oracle, oracle.param.disjunctive_cut_append_rule)
        @test haskey(oracle.dcglp, :con_disjunctive)
        @test length(oracle.dcglp[:con_disjunctive]) == 2

        oracle_none = make_polar_oracle(data, master; append_rule = NoDisjunctiveCuts())
        PolarDCGLP.add_disjunctive_cuts!(oracle_none, oracle_none.param.disjunctive_cut_append_rule)
        @test !haskey(oracle_none.dcglp, :con_disjunctive)

        oracle_all = make_polar_oracle(data, master; append_rule = AllDisjunctiveCuts())
        push!(oracle_all.splits, (sparsevec([1], [1.0], master.dim_x), 0.0))
        PolarDCGLP.append_current_disjunctive_cut!(oracle_all, cut_1)
        PolarDCGLP.add_disjunctive_cuts!(oracle_all, oracle_all.param.disjunctive_cut_append_rule)
        @test haskey(oracle_all.dcglp, :con_disjunctive)
    end

    @testset "Direct Cut Generation" begin
        master = Master(data; customize = customize_master_polar!)
        oracle = make_polar_oracle(data, master)
        is_in_L, _, _ = BendersX.generate_cuts(oracle, fill(0.5, master.dim_x), [0.0]; time_limit = 10.0)
        @test !is_in_L
        @test !isempty(oracle.disjunctiveCuts)
        cut = last(oracle.disjunctiveCuts)
        @test length(cut.a_t) == 1
        @test isapprox(cut.a_t[1], -1.0; atol = 1e-9)
    end

    @testset "BendersSeq Solve Path" begin
        for reuse_dcglp in (true, false), adjust_t_to_fx in (true, false)
            master = Master(data; customize = customize_master_polar!)
            oracle = make_polar_oracle(data, master; reuse_dcglp = reuse_dcglp, adjust_t_to_fx = adjust_t_to_fx)
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

    @testset "Integration Smoke Tests" begin
        master = Master(data; customize = customize_master_polar!)
        typical = ClassicalOracle(data, master; customize = customize_sub_polar!)
        polar = make_polar_oracle(data, master)

        callback = UserCallback(polar)
        @test callback.oracle === polar

        preprocessing = DisjunctiveRootNodePreprocessing(
            typical,
            polar;
            params = BendersSeqParam(time_limit = 10.0, gap_tolerance = 1e-4, verbose = false),
        )
        @test BendersX.root_node_processing!(master, preprocessing) >= 0.0
    end
end
