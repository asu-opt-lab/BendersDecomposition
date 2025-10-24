module BendersBase

# Import dependencies
using JuMP
using MathOptInterface
using Printf
using LinearAlgebra
using SparseArrays
const MOI = MathOptInterface
using DataFrames

# Include source files
include("types.jl")
include("utils/utils.jl")
include("modules/modules.jl") 
include("algorithms/algorithms.jl")


# Export abstract types
export AbstractBendersDecomposition
export Data, AbstractData
export AbstractMaster, AbstractMip
export AbstractOracle, AbstractOracleParam, AbstractTypicalOracle, BasicOracleParam
export AbstractBendersSeq, AbstractBendersCallback
export Seq, SeqInOut
export TerminationStatus, NotSolved, TimeLimit, Optimal, InfeasibleOrNumericalIssue
export TimeLimitException, UnexpectedModelStatusException, UndefError, AlgorithmException
export AbstractLoopState, AbstractLoopLog, AbstractLoopParam
export AbstractBendersSeqState, AbstractBendersSeqLog, AbstractBendersSeqParam
export AbstractBendersBnBState, AbstractBendersBnBLog, AbstractBendersBnBParam
export BendersSeqState, BendersSeqLog, BendersSeqParam
export BendersBnBState, BendersBnBLog, BendersBnBParam
export Hyperplane, aggregate, generate_cuts, set_parameter!

end # module BendersBase