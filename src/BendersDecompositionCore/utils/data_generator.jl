export generate_capacited_facility_location, write_capacited_facility_location_problem, generate_stochastic_capacited_facility_location, write_stochastic_capacited_facility_location_problem

function generate_capacited_facility_location(
    n_facilities::Int,
    n_customers::Int,
    ratio::Int
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
    open("data/CFLP/random_data/$(filename)", "w") do f
        write(f, json_data)
    end
end

# ============================================================================
# SCFLP
# ============================================================================

function generate_stochastic_capacited_facility_location(
    n_facilities::Int,
    n_customers::Int,
    n_scenarios::Int,
    ratio::Int
)

    c_x = rand(n_customers)
    c_y = rand(n_customers)

    f_x = rand(n_facilities)
    f_y = rand(n_facilities)
    

    base_demands = rand(5:35, n_customers)
    
    demand_stds = rand(n_customers) .* (0.02 * base_demands) .+ (0.01 * base_demands)
    

    demands = Vector{Vector{Int}}(undef, n_scenarios)
    # for s in 1:n_scenarios
    #     demands[s] = base_demands
    # end
    demands[1] = base_demands
    for s in 2:n_scenarios
        demands[s] = zeros(Int, n_customers)
        for j in 1:n_customers
            demands[s][j] = max(1, round(Int, rand(Normal(base_demands[j], demand_stds[j]))))
        end #make sure demands are positive
    end

    capacities = rand(10:160, n_facilities)
    fixed_costs = (rand(100:110, n_facilities) .* sqrt.(capacities)) .+ rand(0:90, n_facilities)
    fixed_costs = round.(Int, fixed_costs)


    total_demand_max = maximum([sum(d) for d in demands])
    # total_demand_max = sum(base_demands)
    total_capacity = sum(capacities)

    # adjust capacities according to ratio
    capacities = capacities .* ratio .* total_demand_max ./ total_capacity
    capacities = round.(Int, capacities)


    trans_costs = sqrt.((c_x .- f_x') .^ 2 .+ (c_y .- f_y') .^ 2) .* 10 .* base_demands

    return SCFLPData(n_facilities, n_customers, n_scenarios, capacities, demands, fixed_costs, trans_costs)

end

function write_stochastic_capacited_facility_location_problem(data::SCFLPData; filename::String="data.json")
    json_data = JSON.json(data)
    open("data/SCFLP/$(filename)", "w") do f
        write(f, json_data)
    end
end