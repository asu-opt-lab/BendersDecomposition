using JSON
import BendersX: AbstractData

"""
    SCFLPBodurData <: AbstractData

Data container for the stochastic capacitated facility location problem
using the Bodur-generated instance data. The model methods for this data type
use the same fraction-based second-stage variables as `BendersX.SCFLPData`.
"""
struct SCFLPBodurData <: AbstractData
    n_facilities::Int
    n_customers::Int
    n_scenarios::Int
    capacities::Vector{Float64}
    demands::Vector{Vector{Float64}}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

"""
    read_scflp_bodur(filename; filepath) -> SCFLPBodurData

Read a Bodur-format SCFLP instance from JSON. The JSON schema is identical to
the one used by `read_stochastic_capacited_facility_location_problem`.
"""
function read_scflp_bodur(filename::String; filepath::AbstractString)
    fullpath = joinpath(filepath, "$(filename).json")
    data = JSON.parse(read(fullpath, String))
    costs = reduce(hcat, data["costs"])'
    return SCFLPBodurData(
        data["n_facilities"], data["n_customers"], data["n_scenarios"],
        data["capacities"], data["demands"], data["fixed_costs"], costs,
    )
end
