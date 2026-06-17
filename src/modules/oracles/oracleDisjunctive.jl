"""
Abstract supertype for split disjunctive oracles implemented with a DCGLP model.
"""
abstract type AbstractSplitOracle <: AbstractDisjunctiveOracle end

"""
Prototype for `generate_cuts` on disjunctive oracles.

Concrete subtypes of `AbstractDisjunctiveOracle` should implement this
method to generate disjunctive cuts or fall back to problem-specific typical
cuts when appropriate.
"""
function generate_cuts(oracle::AbstractDisjunctiveOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; time_limit = 3600)
    throw(UndefError("update generate_cuts for $(typeof(oracle))"))
end

function set_parameter!(oracle::AbstractDisjunctiveOracle, args...)
    throw(ArgumentError(
        "set_parameter! is not permitted for AbstractDisjunctiveOracle because their " *
        "parameters must be fixed at construction. Please supply all parameters " *
        "when creating the disjunctive oracle."
    ))
end

include("normalization/normalization.jl")
include("oracleDisjunctiveSplit.jl")
include("oracleDisjunctiveUtils.jl")
include("normalization/normalizationLp.jl")
include("normalization/normalizationReversePolar.jl")
