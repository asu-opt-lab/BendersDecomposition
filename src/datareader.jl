export read_GK_data, read_cflp_benchmark_data, read_uflp_benchmark_data, read_Simple_data, read_Orlib_data


function read_GK_data(filename::AbstractString;filepath="data/random_data/"::AbstractString)
    fullpath = joinpath(filepath, join([filename, ".json"]))
    loaded_json = open(fullpath, "r") do file
       read(file, String)
    end
    loaded_data_string = JSON.parse(loaded_json)
    n_facilities = loaded_data_string["n_facilities"]
    n_customers = loaded_data_string["n_customers"]
    capacities = loaded_data_string["capacities"]
    demands = loaded_data_string["demands"]
    fixed_costs = loaded_data_string["fixed_costs"]
    costs = loaded_data_string["costs"]
    costs = reduce(hcat,costs)'
    return CFLPData(n_facilities, n_customers, capacities, demands, fixed_costs, costs)
end

function read_cflp_benchmark_data(filename::AbstractString;filepath="data/locssall/"::AbstractString)
    fullpath = joinpath(filepath, filename)
    f = open(fullpath)

    line1 = readline(f)
    vals1 = split(line1)
    n_facilities = parse(Int, vals1[1])
    n_customers = parse(Int, vals1[2])

    capacities = zeros(Float32,n_facilities)
    fixed_costs = zeros(Float32,n_facilities)
    for i in 1:n_facilities
        line = readline(f)
        vals = split(line)
        capacities[i] = parse(Float32, vals[1])
        fixed_costs[i] = parse(Float32, vals[2])
    end

    demands = zeros(Float32,n_customers)
    for i in 1:Int(n_customers/10)
        line = readline(f)
        vals = split(line)
        for j in 1:10
            demands[10*(i-1)+j] = parse(Float32, vals[j])
        end
    end

    costs = zeros(Float32,n_facilities,n_customers)
    line_facility = Int(n_customers/10)


    nline = 0
    fth = 1
    while !eof(f)
        
        line = readline(f)
        vals = split(line)
        nline += 1
        for j in 1:10
            costs[fth,10*(nline-1)+j] = parse(Float32, vals[j])
        end
        
        if nline == line_facility
            fth += 1
            nline = 0
        end
        
    end
    
    return CFLPData(n_facilities, n_customers, capacities, demands, fixed_costs, costs)
end

function read_uflp_benchmark_data(filename::AbstractString;filepath="data/locssall/"::AbstractString)
    fullpath = joinpath(filepath, filename)
    f = open(fullpath)

    line1 = readline(f)
    vals1 = split(line1)
    n_facilities = parse(Int, vals1[1])
    n_customers = parse(Int, vals1[2])

    capacities = zeros(Float32,n_facilities)
    fixed_costs = zeros(Float32,n_facilities)
    for i in 1:n_facilities
        line = readline(f)
        vals = split(line)
        capacities[i] = parse(Float32, vals[1])
        fixed_costs[i] = parse(Float32, vals[2])
    end

    demands = zeros(Float32,n_customers)
    for i in 1:Int(n_customers/10)
        line = readline(f)
        vals = split(line)
        for j in 1:10
            demands[10*(i-1)+j] = parse(Float32, vals[j])
        end
    end

    costs = zeros(Float32,n_facilities,n_customers)
    line_facility = Int(n_customers/10)


    nline = 0
    fth = 1
    while !eof(f)
        
        line = readline(f)
        vals = split(line)
        nline += 1
        for j in 1:10
            costs[fth,10*(nline-1)+j] = parse(Float32, vals[j])
        end
        
        if nline == line_facility
            fth += 1
            nline = 0
        end
        
    end
    
    return UFLPData(n_facilities, n_customers, demands, fixed_costs, costs)
end

function read_Simple_data(filename::AbstractString;filepath="data"::AbstractString)
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

function read_Orlib_data(filename::String;filepath="data/M/R"::AbstractString)
    fullpath = joinpath(filepath, filename)
    open(fullpath, "r") do file
        # Read the first line to get n and m
        line = readline(file)
        n, m = parse.(Int, split(line))

        # Initialize arrays to store capacities, opening costs, demands, and allocation costs
        capacities = zeros(Int, n)
        opening_costs = zeros(Float64, n)
        demands = zeros(Int, m)
        allocation_costs = zeros(Float64, m, n)

        # Read the capacities and opening costs for each facility
        for i in 1:n
            line = readline(file)
            capacity, opening_cost = parse.(Float64, split(line))
            capacities[i] = capacity
            opening_costs[i] = opening_cost
        end

        # Read the demands and allocation costs for each city
        for j in 1:m
            demands[j] = parse(Int, readline(file))
            allocation_costs[j, :] = parse.(Float64, split(readline(file)))
        end
        
        if capacities[1] != 0
            return CFLPData(n, m, capacities, demands, opening_costs, allocation_costs)
        else
            return UFLPData(n, m, demands, opening_costs, allocation_costs)
        end
    end
end