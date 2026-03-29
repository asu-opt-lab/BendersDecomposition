module BendersXGurobiExt

using BendersX
using Gurobi
using JuMP

const GRB_ENV = Ref{Union{Nothing,Gurobi.Env}}(nothing)

function shared_gurobi_env()
    if isnothing(GRB_ENV[])
        GRB_ENV[] = Gurobi.Env()
    end
    return GRB_ENV[]::Gurobi.Env
end

function _gurobi_optimizer(args::Pair...)
    return JuMP.optimizer_with_attributes(
        () -> Gurobi.Optimizer(shared_gurobi_env()),
        BendersX._normalize_optimizer_attributes(args...)...,
    )
end

function __init__()
    BendersX._register_gurobi_optimizer!(_gurobi_optimizer)
    return nothing
end

end
