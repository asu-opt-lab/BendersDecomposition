"""
Abstract type for disjunctive oracles
"""
abstract type AbstractDisjunctiveOracle <: AbstractOracle end

function generate_cuts(oracle::AbstractDisjunctiveOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; time_limit = 3600)
    throw(UndefError("update generate_cuts for $(typeof(oracle))"))
end

# Common utility functions for managing oracle parameters
function set_parameter!(oracle::AbstractDisjunctiveOracle, args...)
    throw(ArgumentError(
        "set_parameter! is not permitted for DisjunctiveOracles because their " *
        "parameters must be fixed at construction. Please supply all parameters " *
        "when creating the disjunctive oracle."
    ))
end

include("oracleDisjunctive.jl")
