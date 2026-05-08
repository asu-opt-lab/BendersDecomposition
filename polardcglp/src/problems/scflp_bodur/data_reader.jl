using JSON
import BendersX: AbstractData

"""
    SCFLPBodurData <: AbstractData

Data container for the stochastic capacitated facility location problem
following the Bodur et al. (2016) formulation (eq. 15a-e), where second-stage
variables `y[i,j,k]` represent absolute flow (amount of customer j's demand
served by facility i in scenario k), as opposed to the fraction-based
formulation in `BendersX.SCFLPData`.
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
