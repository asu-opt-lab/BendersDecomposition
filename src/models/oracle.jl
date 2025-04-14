export AbstractTypicalOracle, AbstractDisjunctiveOracle, generate_cuts

abstract type AbstractTypicalOracle <: AbstractOracle end
abstract type AbstractDisjunctiveOracle <: AbstractOracle end

function generate_cuts(oracle::AbstractOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol = 1e-6, time_limit = 3600)
    throw(UndefError("update generate_cuts for $(typeof(AbstractOracle))"))
end

include("oracleTypical.jl")
include("oracleTypicalSeparable.jl")
include("oracleDisjunctive.jl")