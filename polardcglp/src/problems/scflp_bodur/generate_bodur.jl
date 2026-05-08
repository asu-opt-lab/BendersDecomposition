using Random
using JSON

# Reuses read_flcap_data and write_scflp_json from the PolarDCGLP problem includes.

"""
    generate_scflp_bodur(flcap_filename, n_scenarios; seed, sigma_l, sigma_u, flcap_dir)
        -> Dict

Generate a stochastic SCFLP instance from an FLCAP base file using the strict
Bodur et al. (2016) sampling rule:
- For each customer j, draw σ_j ~ U(sigma_l·μ̄_j, sigma_u·μ̄_j) ONCE.
- For each scenario k, draw ξ_{k,j} ~ N(μ_j, σ_j²), truncated at 0.

Default seed = n_facilities*10000 + n_customers*1000 + n_scenarios.
"""
function generate_scflp_bodur(flcap_filename::AbstractString, n_scenarios::Int;
                              seed::Union{Int,Nothing}=nothing,
                              sigma_l::Float64=0.1,
                              sigma_u::Float64=0.3,
                              flcap_dir::AbstractString=joinpath(POLARDCGLP_ROOT, "data", "FLCAP"))
    cflp_data = read_flcap_data(flcap_filename; filepath=flcap_dir)

    actual_seed = seed === nothing ?
        cflp_data.n_facilities * 10000 + cflp_data.n_customers * 1000 + n_scenarios :
        seed

    rng = MersenneTwister(actual_seed)

    sigma_per_j = Vector{Float64}(undef, cflp_data.n_customers)
    for j in 1:cflp_data.n_customers
        sigma_per_j[j] = sigma_l + (sigma_u - sigma_l) * rand(rng)
    end

    demands_out = Vector{Vector{Float64}}(undef, n_scenarios)
    for k in 1:n_scenarios
        xi_k = Vector{Float64}(undef, cflp_data.n_customers)
        for j in 1:cflp_data.n_customers
            Z = randn(rng)
            xi_k[j] = max(0.0,
                cflp_data.demands[j] + sigma_per_j[j] * cflp_data.demands[j] * Z)
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
