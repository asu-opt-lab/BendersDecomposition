module PolarDCGLP

using JuMP
using MathOptInterface
using LinearAlgebra
using Printf
using Random
using SparseArrays

import BendersX

const MOI = MathOptInterface
const POLAR_T_LOWER_BOUND = -1e6

export PolarDCGLPParam, PolarDCGLPOracle

mutable struct PolarDCGLPParam <: BendersX.AbstractOracleParam
    dcglp_param::BendersX.DcglpParam
    split_index_selection_rule::BendersX.SplitIndexSelectionRule
    disjunctive_cut_append_rule::BendersX.DisjunctiveCutsAppendRule
    add_benders_cuts_to_master::Int
    fraction_of_benders_cuts_to_master::Float64
    reuse_dcglp::Bool
    strengthened::Bool
    lift::Bool
    zero_tol::Float64

    function PolarDCGLPParam(
        dcglp_param::BendersX.DcglpParam;
        split_index_selection_rule::BendersX.SplitIndexSelectionRule = BendersX.RandomFractional(),
        disjunctive_cut_append_rule::BendersX.DisjunctiveCutsAppendRule = BendersX.AllDisjunctiveCuts(),
        add_benders_cuts_to_master::Union{Bool, Int} = 1,
        fraction_of_benders_cuts_to_master::Float64 = 1.0,
        reuse_dcglp::Bool = true,
        strengthened::Bool = true,
        lift::Bool = false,
        zero_tol::Float64 = 1e-9,
    )
        add_bcuts_to_master =
            add_benders_cuts_to_master === true ? 1 :
            add_benders_cuts_to_master === false ? 0 :
            add_benders_cuts_to_master in (0, 1, 2) ? add_benders_cuts_to_master :
            throw(ArgumentError("`add_benders_cuts_to_master` must be true, false, or an integer in {0, 1, 2}."))

        0.0 < fraction_of_benders_cuts_to_master <= 1.0 ||
            throw(ArgumentError("`fraction_of_benders_cuts_to_master` must lie in (0, 1]."))

        new(
            dcglp_param,
            split_index_selection_rule,
            disjunctive_cut_append_rule,
            add_bcuts_to_master,
            fraction_of_benders_cuts_to_master,
            reuse_dcglp,
            strengthened,
            lift,
            zero_tol,
        )
    end
end

mutable struct PolarDCGLPOracle <: BendersX.AbstractDisjunctiveOracle
    param::PolarDCGLPParam
    dcglp::Model
    typical_oracles::Vector{BendersX.AbstractTypicalOracle}
    disjunctiveCutsByIndex::Vector{Vector{BendersX.Hyperplane}}
    disjunctiveCuts::Vector{BendersX.Hyperplane}
    splits::Vector{Tuple{SparseVector{Float64, Int}, Float64}}

    function PolarDCGLPOracle(
        master::BendersX.AbstractMaster,
        typical_oracles::Vector{T},
        param::PolarDCGLPParam,
    ) where {T <: BendersX.AbstractTypicalOracle}
        length(typical_oracles) == 2 ||
            throw(ArgumentError("PolarDCGLPOracle requires exactly two typical oracles."))

        for xi in master.x
            is_binary(xi) || throw(ArgumentError("PolarDCGLPOracle requires all master variables to be binary."))
        end

        dcglp = build_polar_dcglp(master, param)
        disjunctive_cuts_by_index = [Vector{BendersX.Hyperplane}() for _ in 1:master.dim_x]
        splits = Vector{Tuple{SparseVector{Float64, Int}, Float64}}()

        new(param, dcglp, Vector{BendersX.AbstractTypicalOracle}(typical_oracles), disjunctive_cuts_by_index, BendersX.Hyperplane[], splits)
    end
end

function build_polar_dcglp(master::BendersX.AbstractMaster, param::PolarDCGLPParam)
    dcglp = Model(param.dcglp_param.optimizer)

    @variable(dcglp, tau[1:master.dim_t])
    @variable(dcglp, omega_0[1:2] >= 0)
    @variable(dcglp, omega_x[1:2, 1:master.dim_x])
    @variable(dcglp, omega_t[1:2, 1:master.dim_t])

    @objective(dcglp, Min, sum(tau))

    @constraint(dcglp, [i in 1:2], omega_t[i, :] .>= POLAR_T_LOWER_BOUND .* omega_0[i])
    @constraint(dcglp, coneta[i in 1:2, j in 1:master.dim_x], 0 >= -omega_0[i] + omega_x[i, j])
    @constraint(dcglp, condelta[i in 1:2, j in 1:master.dim_x], 0 >= -omega_x[i, j])

    @constraint(dcglp, con0, omega_0[1] + omega_0[2] == 1.0)
    @constraint(dcglp, conx[j in 1:master.dim_x], omega_x[1, j] + omega_x[2, j] == 0.0)
    @constraint(dcglp, cont[j in 1:master.dim_t], omega_t[1, j] + omega_t[2, j] - tau[j] == 0.0)

    for i in 1:2
        BendersX.transfer_scaled_linear_rows_and_bounds_with_types!(master.model, master.x, dcglp, dcglp[:omega_x][i, :], dcglp[:omega_0][i])
    end

    return dcglp
end

function delete_registered_constraints!(model::Model, sym::Symbol)
    haskey(model, sym) || return
    registered = model[sym]
    if registered isa AbstractArray
        delete.(model, registered)
    else
        delete(model, registered)
    end
    unregister(model, sym)
end

function get_split_index(oracle::PolarDCGLPOracle)
    isa(oracle.param.split_index_selection_rule, BendersX.SimpleSplit) ||
        throw(BendersX.AlgorithmException("get_split_index is only valid for simple split rules."))
    return findfirst(x -> x > 0.5, oracle.splits[end][1])
end

function replace_disjunctive_inequality!(oracle::PolarDCGLPOracle)
    dcglp = oracle.dcglp
    phi, phi_0 = oracle.splits[end]

    delete_registered_constraints!(dcglp, :con_split_kappa)
    delete_registered_constraints!(dcglp, :con_split_nu)

    @constraint(dcglp, con_split_kappa, 0 >= dcglp[:omega_0][1] * (phi_0 + 1.0) - phi' * dcglp[:omega_x][1, :])
    @constraint(dcglp, con_split_nu, 0 >= -dcglp[:omega_0][2] * phi_0 + phi' * dcglp[:omega_x][2, :])

    return nothing
end

include("PolarDCGLPInterface.jl")
include("solvePolarDcglp.jl")

end # module PolarDCGLP
