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

include("disjunctive/normalization/normalization.jl")
include("disjunctive/splitOracle.jl")
include("disjunctive/splitUtils.jl")
include("disjunctive/cutUtils.jl")
include("disjunctive/dcglpUtils.jl")
include("disjunctive/loggingUtils.jl")
include("disjunctive/normalization/lpDistance.jl")
include("disjunctive/normalization/reversePolar.jl")
