# This file is part of the BendersDecomposition.jl library.
# Copyright (c) 2024, Kaiwen Fang
# All rights reserved.

module BendersDecomposition

using Printf, StatsBase, Random, Distributions, LinearAlgebra, Plots, ArgParse, DataFrames, CSV, JSON
using JuMP, CPLEX, Gurobi

# Abstract type for problem data
abstract type AbstractData end
include("datasets.jl")

# Abstract types for different problem components
abstract type AbstractMasterProblem end
abstract type AbstractSubProblem end
abstract type AbstractDCGLP end
abstract type AbstractMILP end

# Abstract types for Benders algorithm and cut generation
abstract type AbstractBendersAlgorithm end
abstract type CutGenerationStrategy end
struct StandardCut <: CutGenerationStrategy end
struct FatKnapsackCut <: CutGenerationStrategy end
struct SlimKnapsackCut <: CutGenerationStrategy end
export StandardCut, FatKnapsackCut, SlimKnapsackCut

# Abstract types for different norms
abstract type AbstractNormType end
struct StandardNorm <: AbstractNormType end
abstract type LNorm <: AbstractNormType end

struct L1Norm <: LNorm end
struct L2Norm <: LNorm end
struct LInfNorm <: LNorm end
export StandardNorm, L1Norm, L2Norm, LInfNorm
# Parameters for Benders algorithm
struct BendersParams
    time_limit::Float64
    gap_tolerance::Float64
end
export BendersParams
# Include model and algorithm implementations
include("models/models.jl")
include("algorithms/algorithms.jl")
include("utilize.jl")

end
