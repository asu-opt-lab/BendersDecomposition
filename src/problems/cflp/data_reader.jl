
using JSON
using LinearAlgebra

"""
    CFLPData <: AbstractData

Data container for the capacitated facility location problem.
"""
struct CFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    capacities::Vector{Float64}
    demands::Vector{Float64}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

"""
    read_GK_data(filename; filepath = get_artifact_path("cflp_random_data")) -> CFLPData

Read a JSON-formatted random CFLP instance from the packaged GK dataset.
"""
function read_GK_data(filename::AbstractString; filepath=get_artifact_path("cflp_random_data"))
    fullpath = joinpath(filepath, "$(filename).json")
    loaded_json = read(fullpath, String)
    data = JSON.parse(loaded_json)
    costs = reduce(hcat, data["costs"])'
    return CFLPData(data["n_facilities"], data["n_customers"], 
                    data["capacities"], data["demands"], data["fixed_costs"], costs)
end

"""
    read_cflp_benchmark_data(filename; filepath = get_artifact_path("cflp_locssall")) -> CFLPData

Read a benchmark CFLP instance from the packaged LOCSSALL-style dataset.
"""
function read_cflp_benchmark_data(filename::AbstractString; filepath=get_artifact_path("cflp_locssall"))
    fullpath = joinpath(filepath, filename)

    n_facilities, n_customers, capacities, fixed_costs, demands, costs = open(fullpath) do f
        vals1 = split(readline(f))
        n_facilities = parse(Int, vals1[1])
        n_customers = parse(Int, vals1[2])

        capacities = zeros(Float64, n_facilities)
        fixed_costs = zeros(Float64, n_facilities)
        for i in 1:n_facilities
            vals = split(readline(f))
            capacities[i] = parse(Float64, vals[1])
            fixed_costs[i] = parse(Float64, vals[2])
        end

        tokens = split(read(f, String))
        expected = n_customers + n_facilities * n_customers
        length(tokens) == expected || error("Expected $(expected) remaining numeric values in $(fullpath), found $(length(tokens))")

        demands = zeros(Float64, n_customers)
        idx = 1
        for j in 1:n_customers
            demands[j] = parse(Float64, tokens[idx])
            idx += 1
        end

        costs = zeros(Float64, n_facilities, n_customers)
        for i in 1:n_facilities
            for j in 1:n_customers
                costs[i, j] = parse(Float64, tokens[idx])
                idx += 1
            end
        end

        n_facilities, n_customers, capacities, fixed_costs, demands, costs
    end

    return CFLPData(n_facilities, n_customers, capacities, demands, fixed_costs, costs)
end

"""
    read_cfl_file(filename; filepath = get_artifact_path("cflp_output")) -> CFLPData

Read a `.cfl` formatted CFLP instance from the packaged output dataset.
"""
function read_cfl_file(filename::AbstractString; filepath=get_artifact_path("cflp_output"))
    fullpath = joinpath(filepath, "$(filename).cfl")
    f = open(fullpath)
    
    line = readline(f)
    startswith(line, "[CFLP-PROBLEMFILE]") || error("File format not recognized")
    readline(f)  # timestamp
    
    line = readline(f)
    n_customers = parse(Int, match(r"#customers:\s*(\d+)", line).captures[1])
    n_facilities = parse(Int, match(r"#depot sites:\s*(\d+)", line).captures[1])
    
    readline(f)  # blank
    readline(f)  # [DEPOTS]
    readline(f)  # headers
    
    capacities = zeros(Float64, n_facilities)
    fixed_costs = zeros(Float64, n_facilities)
    for i in 1:n_facilities
        vals = split(readline(f))
        capacities[i] = parse(Float64, vals[1])
        fixed_costs[i] = parse(Float64, vals[2])
    end
    
    readline(f)  # blank
    readline(f)  # [CUSTOMERS]
    readline(f)  # headers
    
    demands = zeros(Float64, n_customers)
    for i in 1:n_customers
        vals = split(readline(f))
        demands[i] = parse(Float64, vals[1])
    end
    
    readline(f)  # blank
    readline(f)  # [COSTMATRIX]
    readline(f)  # formula
    readline(f)  # [MATRIX]
    readline(f)  # Dim line
    
    costs = zeros(Float64, n_facilities, n_customers)
    for i in 1:n_facilities
        vals = split(readline(f))
        for j in 1:n_customers
            costs[i, j] = parse(Float64, vals[j]) / demands[j]
        end
    end
    close(f)
    
    return CFLPData(n_facilities, n_customers, capacities, demands, fixed_costs, costs)
end
