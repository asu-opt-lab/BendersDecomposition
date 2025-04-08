# export UFLPData, read_uflp_benchmark_data, read_Simple_data

struct UFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    demands::Vector{Float64}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

function read_uflp_benchmark_data(filename::AbstractString;filepath="example/uflp/data/locssall/"::AbstractString)
    fullpath = joinpath(filepath, filename)
    f = open(fullpath)

    line1 = readline(f)
    vals1 = split(line1)
    n_facilities = parse(Int, vals1[1])
    n_customers = parse(Int, vals1[2])

    capacities = zeros(Float64,n_facilities)
    fixed_costs = zeros(Float64,n_facilities)
    for i in 1:n_facilities
        line = readline(f)
        vals = split(line)
        capacities[i] = parse(Float64, vals[1])
        fixed_costs[i] = parse(Float64, vals[2])
    end

    demands = zeros(Float64,n_customers)
    for i in 1:Int(n_customers/10)
        line = readline(f)
        vals = split(line)
        for j in 1:10
            demands[10*(i-1)+j] = parse(Float64, vals[j])
        end
    end

    costs = zeros(Float64,n_facilities,n_customers)
    line_facility = Int(n_customers/10)


    nline = 0
    fth = 1
    while !eof(f)
        
        line = readline(f)
        vals = split(line)
        nline += 1
        for j in 1:10
            costs[fth,10*(nline-1)+j] = parse(Float64, vals[j])
        end
        
        if nline == line_facility
            fth += 1
            nline = 0
        end
        
    end
    
    return UFLPData(n_facilities, n_customers, demands, fixed_costs, costs)
end

function read_Simple_data(filename::AbstractString;filepath="example/uflp/data/AllKoerkelGhosh"::AbstractString)
    fullpath = joinpath(filepath, filename)
    f = open(fullpath)

    readline(f)
    line1 = readline(f)
    vals1 = split(line1)
    n_facilities = parse(Int, vals1[1])
    n_customers = parse(Int, vals1[2])

    fixed_costs = zeros(Int,n_facilities)
    costs = zeros(Int,n_facilities,n_customers)

    fth = 1
    while !eof(f)
        line = readline(f)
        vals = split(line)
        
        fixed_costs[fth] = parse(Int, vals[2])
        for j in 1:n_customers
            costs[fth,j] = parse(Int, vals[2+j])
        end
        fth += 1
    end

    demands = ones(Int,n_customers)
    return UFLPData(n_facilities, n_customers, demands, fixed_costs, costs)
end