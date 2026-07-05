
"""
    UFLPData <: AbstractData

Data container for the uncapacitated facility location problem.
"""
struct UFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    demands::Vector{Float64}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

"""
    read_uflp_benchmark_data(filename; filepath = get_artifact_path("uflp_locssall")) -> UFLPData

Read a benchmark UFLP instance from the packaged LOCSSALL-style dataset.
"""
function read_uflp_benchmark_data(filename::AbstractString; filepath=get_artifact_path("uflp_locssall"))
    fullpath = joinpath(filepath, filename)

    n_facilities, n_customers, demands, fixed_costs, costs = open(fullpath) do f
        vals1 = split(readline(f))
        n_facilities = parse(Int, vals1[1])
        n_customers = parse(Int, vals1[2])

        fixed_costs = zeros(Float64, n_facilities)
        for i in 1:n_facilities
            vals = split(readline(f))
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

        n_facilities, n_customers, demands, fixed_costs, costs
    end

    return UFLPData(n_facilities, n_customers, demands, fixed_costs, costs)
end

"""
    read_Simple_data(filename; filepath = get_artifact_path("uflp_allkoerkelghosh")) -> UFLPData

Read a UFLP instance from the packaged "Simple" dataset.
"""
function read_Simple_data(filename::AbstractString; filepath=get_artifact_path("uflp_allkoerkelghosh"))
    fullpath = joinpath(filepath, filename)
    f = open(fullpath)

    readline(f)
    vals1 = split(readline(f))
    n_facilities = parse(Int, vals1[1])
    n_customers = parse(Int, vals1[2])

    fixed_costs = zeros(Int, n_facilities)
    costs = zeros(Int, n_facilities, n_customers)

    fth = 1
    while !eof(f)
        vals = split(readline(f))
        fixed_costs[fth] = parse(Int, vals[2])
        for j in 1:n_customers
            costs[fth, j] = parse(Int, vals[2+j])
        end
        fth += 1
    end
    close(f)

    return UFLPData(n_facilities, n_customers, ones(Int, n_customers), fixed_costs, costs)
end
