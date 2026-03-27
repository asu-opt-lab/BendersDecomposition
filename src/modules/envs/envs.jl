include("BendersSeq.jl")
include("BendersBnB.jl")
include("BendersSeqInOut.jl")
include("SpecializedBendersSeq.jl")
include("callback/preprocessingDisjunctive.jl")

"""
    solve!(env::AbstractBendersEnv)

Run a configured Benders environment and return its execution log.

Concrete environment types such as [`BendersSeq`](@ref), [`BendersSeqInOut`](@ref),
[`BendersBnB`](@ref), and [`SpecializedBendersSeq`](@ref) provide the actual
implementations. This fallback exists to document the public interface and to
raise a descriptive error for unsupported environment types.
"""
function solve!(env::AbstractBendersEnv)
    throw(UndefError("update solve! for $(typeof(env))"))
end
