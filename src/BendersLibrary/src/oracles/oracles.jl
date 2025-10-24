"""
Abstract type for disjunctive oracles
"""
abstract type AbstractDisjunctiveOracle <: AbstractOracle end

function generate_cuts(oracle::AbstractDisjunctiveOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; time_limit = 3600)
    throw(UndefError("update generate_cuts for $(typeof(oracle))"))
end

include("oracleDisjunctive.jl")
