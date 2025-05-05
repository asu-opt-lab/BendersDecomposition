# export read_GK_data, read_cflp_benchmark_data

using JSON
using LinearAlgebra
struct CFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    capacities::Vector{Float64}
    demands::Vector{Float64}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

function read_GK_data(filename::AbstractString;filepath="example/cflp/data/random_data/"::AbstractString)
    fullpath = joinpath(filepath, join([filename, ".json"]))
    # fullpath = joinpath(filepath, filename)
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

function read_cflp_benchmark_data(filename::AbstractString;filepath="example/cflp/data/locssall/"::AbstractString)
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
    
    return CFLPData(n_facilities, n_customers, capacities, demands, fixed_costs, costs)
end

function read_cfl_file(filename::AbstractString; filepath="example/cflp/data/output/"::AbstractString)
    fullpath = joinpath(filepath, join([filename, ".cfl"]))
    f = open(fullpath)
    
    # Skip the header line [CFLP-PROBLEMFILE]
    line = readline(f)
    if !startswith(line, "[CFLP-PROBLEMFILE]")
        close(f)
        error("File format not recognized as CFLP format")
    end
    
    # Skip the timestamp line
    readline(f)
    
    # Read problem dimensions and ratio
    line = readline(f)
    m_str = match(r"#customers:\s*(\d+)", line)
    n_str = match(r"#depot sites:\s*(\d+)", line)
    
    n_customers = parse(Int, m_str.captures[1])
    n_facilities = parse(Int, n_str.captures[1])
    
    # Skip blank line
    readline(f)
    
    # Skip [DEPOTS] header
    line = readline(f)
    if !startswith(line, "[DEPOTS]")
        close(f)
        error("Expected [DEPOTS] section not found")
    end
    
    # Skip column headers
    readline(f)
    
    # Read facility data
    capacities = zeros(Float64, n_facilities)
    fixed_costs = zeros(Float64, n_facilities)
    facility_coords = zeros(Float64, n_facilities, 2)
    
    for i in 1:n_facilities
        line = readline(f)
        vals = split(line)
        @assert length(vals) == 6
        capacities[i] = parse(Float64, vals[1])
        fixed_costs[i] = parse(Float64, vals[2])
        facility_coords[i, 1] = parse(Float64, vals[4])
        facility_coords[i, 2] = parse(Float64, vals[5])
    end
    
    # Skip blank line
    readline(f)
    
    # Skip [CUSTOMERS] header
    line = readline(f)
    if !startswith(line, "[CUSTOMERS]")
        close(f)
        error("Expected [CUSTOMERS] section not found")
    end
    
    # Skip column headers
    readline(f)
    
    # Read customer data
    demands = zeros(Float64, n_customers)
    customer_coords = zeros(Float64, n_customers, 2)
    
    for i in 1:n_customers
        line = readline(f)
        vals = split(line)
        @assert length(vals) == 4
        demands[i] = parse(Float64, vals[1])
        customer_coords[i, 1] = parse(Float64, vals[2])
        customer_coords[i, 2] = parse(Float64, vals[3])
    end
    
    # Skip blank line
    readline(f)
    
    # Skip [COSTMATRIX] header and formula line
    line = readline(f)
    if !startswith(line, "[COSTMATRIX]")
        close(f)
        error("Expected [COSTMATRIX] section not found")
    end
    readline(f)
    
    # Skip [MATRIX] header
    line = readline(f)
    if !startswith(line, "[MATRIX]")
        close(f)
        error("Expected [MATRIX] section not found")
    end
    
    # Read dimension line and verify dimensions
    line = readline(f)
    dim_match = match(r"Dim\s+(\d+)\s+(\d+)", line)
    if dim_match === nothing || 
       parse(Int, dim_match.captures[1]) != n_facilities || 
       parse(Int, dim_match.captures[2]) != n_customers
        close(f)
        error("Matrix dimensions in file do not match expected dimensions")
    end
    
    # Read cost matrix
    costs = zeros(Float64, n_facilities, n_customers)
    
    for i in 1:n_facilities
        line = readline(f)
        vals = split(line)
        @assert length(vals) == n_customers
        for j in 1:n_customers
            costs[i, j] = parse(Float64, vals[j]) #/ demands[j]
        end
    end
    
    # costs2 = zeros(Float64, n_customers, n_facilities)
    # for i in 1:n_facilities
    #     for j in 1:n_customers
    #         costs2[i, j] = round(norm(facility_coords[i, :] - customer_coords[j, :],2) * 0.01 * demands[j], digits=4)
    #     end
    #     println("facility $i: ", costs2[i, :])
    # end
    

    close(f)
    
    return CFLPData(n_facilities, n_customers, capacities, demands, fixed_costs, costs)
end