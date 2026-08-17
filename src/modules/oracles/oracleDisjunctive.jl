"""
Abstract supertype for oracles generating split cuts.
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

include("disjunctive/utilsLoopDcglp.jl")
include("disjunctive/normalization/normalization.jl")
include("disjunctive/Oraclesplit.jl")
include("disjunctive/utilsSplit.jl")
include("disjunctive/utilsCuts.jl")
include("disjunctive/utilsDcglp.jl")
include("disjunctive/utilsLogging.jl")
include("disjunctive/normalization/lpDistance.jl")
include("disjunctive/normalization/reversePolar.jl")
