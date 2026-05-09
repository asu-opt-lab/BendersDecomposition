module SimplexNormDCGLP

using JuMP
using MathOptInterface
using LinearAlgebra
using Printf
using Random
using SparseArrays
using JSON

import BendersX
import BendersX: customize_mip_model!, customize_master_model!, customize_sub_model!

const MOI = MathOptInterface
const POLARDCGLP_ROOT = normpath(joinpath(@__DIR__, ".."))
const SIMPLEXNORM_T_LOWER_BOUND = -1e6
const VERTICAL_REVERSE_POLAR_T_LOWER_BOUND = -1e6
const SIMPLEXNORMTEST_T_LOWER_BOUND = -1e6

export SimplexNormDCGLPParam, SimplexNormDCGLPOracle
export DirectionalPolarDCGLPParam, DirectionalPolarDCGLPOracle, set_core_point!
export VerticalReversePolarDCGLPParam, VerticalReversePolarDCGLPOracle
export SimplexNormTestDCGLPParam, SimplexNormTestDCGLPOracle
export SCFLPBodurData, read_flcap_data, read_scflp_bodur
export generate_scflp_from_flcap, generate_scflp_bodur, write_scflp_json

include("modules/modules.jl")
include("problems/problems.jl")

end # module SimplexNormDCGLP
