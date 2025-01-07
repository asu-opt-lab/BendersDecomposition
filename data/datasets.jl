include("DataReader.jl")

abstract type AbstractData end

struct CFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    capacities::Vector{Float64}
    demands::Vector{Float64}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

struct UFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    demands::Vector{Float64}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

function generate_capacited_facility_location(
    n_facilities::Int,
    n_customers::Int,
    ratio::Float64
)

    c_x = rand(n_customers)
    c_y = rand(n_customers)

    f_x = rand(n_facilities)
    f_y = rand(n_facilities)

    demands = rand(5:35, n_customers)
    capacities = rand(10:160, n_facilities)
    fixed_costs = (rand(100:110, n_facilities) .* sqrt.(capacities)) .+ rand(0:90, n_facilities)
    fixed_costs = round.(Int, fixed_costs)

    total_demand = sum(demands)
    total_capacity = sum(capacities)

    # adjust capacities according to ratio
    capacities = capacities .* ratio .* total_demand ./ total_capacity
    capacities = round.(Int, capacities)
    total_capacity = sum(capacities)

    # transportation costs
    trans_costs = sqrt.((c_x .- f_x') .^ 2 .+ (c_y .- f_y') .^ 2) .* 10 .* demands
    

    return CFLPData(n_facilities, n_customers, capacities, demands, fixed_costs, trans_costs)

end


function write_capacited_facility_location_problem(data::CFLPData; filename::String="data.json")
    json_data = JSON.json(data)
    open("src/BendersDatasets/random_data/$(filename)", "w") do f
        write(f, json_data)
    end
end