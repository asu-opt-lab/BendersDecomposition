abstract type AbstractBendersSeq <: AbstractBendersDecomposition end
abstract type AbstractBendersCallback <: AbstractBendersDecomposition end

include("BendersSeq.jl") 
include("BendersSeqInOut.jl") 
include("BendersBnB.jl") # must be after BendersSeq files
include("Dcglp.jl") 
include("SpecializedBendersSeq.jl") 
