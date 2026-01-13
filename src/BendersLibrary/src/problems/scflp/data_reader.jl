export SCFLPData, read_stochastic_capacited_facility_location_problem
using JSON

struct SCFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    n_scenarios::Int
    capacities::Vector{Float64}
    demands::Vector{Vector{Float64}}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

function read_stochastic_capacited_facility_location_problem(filename::String; filepath=get_artifact_path("scflp"))
    fullpath = joinpath(filepath, "$(filename).json")
    data = JSON.parse(read(fullpath, String))
    costs = reduce(hcat, data["costs"])'
    return SCFLPData(data["n_facilities"], data["n_customers"], data["n_scenarios"],
                     data["capacities"], data["demands"], data["fixed_costs"], costs)
end
