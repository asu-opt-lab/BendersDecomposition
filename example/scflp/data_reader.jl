using JSON

struct SCFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    n_scenarios::Int
    capacities::Vector{Float64}
    demands::Vector{Vector{Float64}} # multiple scenarios
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

function read_stochastic_capacited_facility_location_problem(filename::String;filepath="example/scflp/data/SCFLP/"::AbstractString)
    fullpath = joinpath(filepath, join([filename, ".json"]))
    loaded_json = open(fullpath, "r") do file
        read(file, String)
    end
    loaded_data_string = JSON.parse(loaded_json)
    n_facilities = loaded_data_string["n_facilities"]
    n_customers = loaded_data_string["n_customers"]
    n_scenarios = loaded_data_string["n_scenarios"]
    capacities = loaded_data_string["capacities"]
    demands = loaded_data_string["demands"]
    fixed_costs = loaded_data_string["fixed_costs"]
    costs = loaded_data_string["costs"]
    costs = reduce(hcat,costs)'
    return SCFLPData(n_facilities, n_customers, n_scenarios, capacities, demands, fixed_costs, costs)
end

