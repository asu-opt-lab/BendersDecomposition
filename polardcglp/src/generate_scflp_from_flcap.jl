using Random
using JSON

function generate_scflp_from_flcap(flcap_filename::AbstractString, n_scenarios::Int;
                                    seed::Union{Int,Nothing}=nothing,
                                    sigma_l::Float64=0.1,
                                    sigma_u::Float64=0.3,
                                    flcap_dir::AbstractString=joinpath(@__DIR__, "..", "data", "FLCAP"))
    cflp_data = read_flcap_data(flcap_filename; filepath=flcap_dir)

    actual_seed = seed === nothing ?
        cflp_data.n_customers * 10000 + cflp_data.n_facilities * 1000 + n_scenarios :
        seed

    rng = MersenneTwister(actual_seed)

    demands_out = Vector{Vector{Float64}}(undef, n_scenarios)

    for k in 1:n_scenarios
        xi_k = Vector{Float64}(undef, cflp_data.n_customers)
        for j in 1:cflp_data.n_customers
            sigma = sigma_l + (sigma_u - sigma_l) * rand(rng)
            Z = randn(rng)
            xi_k[j] = max(0.0, cflp_data.demands[j] * (1.0 + sigma * Z))
        end
        demands_out[k] = xi_k
    end

    costs_rows = Vector{Vector{Float64}}(undef, cflp_data.n_facilities)
    for i in 1:cflp_data.n_facilities
        costs_rows[i] = collect(cflp_data.costs[i, :])
    end

    return Dict(
        "n_facilities" => cflp_data.n_facilities,
        "n_customers" => cflp_data.n_customers,
        "n_scenarios" => n_scenarios,
        "capacities" => collect(cflp_data.capacities),
        "demands" => demands_out,
        "fixed_costs" => collect(cflp_data.fixed_costs),
        "costs" => costs_rows,
    )
end

function write_scflp_json(payload::Dict, output_path::AbstractString)
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        JSON.print(io, payload)
    end
    return nothing
end
