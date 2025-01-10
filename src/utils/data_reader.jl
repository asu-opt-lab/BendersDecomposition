export read_GK_data, read_cflp_benchmark_data, read_uflp_benchmark_data, read_Simple_data, read_Orlib_data, read_stochastic_capacited_facility_location_problem
export read_mcndp_instance, read_snip_data

function read_GK_data(filename::AbstractString;filepath="data/CFLP/random_data/"::AbstractString)
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


function read_stochastic_capacited_facility_location_problem(filename::String;filepath="data/SCFLP/"::AbstractString)
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


function read_mcndp_instance(filename::String;filepath="data/NDR/"::AbstractString)
    fullpath = joinpath(filepath, filename)
    open(fullpath, "r") do f
        # Skip the filename line
        readline(f)
        
        # Read problem dimensions
        dims = split(readline(f))
        num_nodes = parse(Int, dims[1])
        num_arcs = parse(Int, dims[2])
        num_commodities = parse(Int, dims[3])
        
        # Initialize data structures
        arcs = Tuple{Int,Int}[]
        fixed_costs = Float64[]
        variable_costs = Float64[]
        capacities = Float64[]
        demands = Tuple{Int,Int,Float64}[]
        
        # Read main data section: arc information
        for i in 1:num_arcs
            line = split(readline(f))
            i_from = parse(Int, line[1])
            i_to = parse(Int, line[2])
            fixed = parse(Float64, line[5])
            var_cost = parse(Float64, line[3])
            capacity = parse(Float64, line[4])
            
            push!(arcs, (i_from, i_to))
            push!(fixed_costs, fixed)
            push!(variable_costs, var_cost)
            push!(capacities, capacity)
        end
        
        # Read commodity demand information
        while !eof(f)
            line = split(readline(f))
            if length(line) >= 3
                origin = parse(Int, line[1])
                dest = parse(Int, line[2])
                demand = parse(Float64, line[3])
                push!(demands, (origin, dest, demand))
            end
        end
        
        return MCNDPData(
            num_nodes,
            num_arcs,
            num_commodities,
            arcs,
            fixed_costs,
            variable_costs,
            capacities,
            demands
        )
    end
end

function create_node_mapping(D, A_minus_D, S)
    # Collect all nodes that appear in the network
    nodes = Set{Int}()
    
    # Collect nodes from D (sensor installation arcs)
    for (i, j, _, _) in D
        push!(nodes, i, j)
    end
    
    # Collect nodes from A_minus_D (non-sensor arcs)
    for (i, j, _) in A_minus_D
        push!(nodes, i, j)
    end
    
    # Collect nodes from scenarios S
    for (i, j, _) in S
        push!(nodes, i, j)
    end
    
    # Create mapping from old to new node indices
    sorted_nodes = sort(collect(nodes))
    node_mapping = Dict(old => new for (new, old) in enumerate(sorted_nodes))
    
    return node_mapping
end

function read_snip_data(instance_no::Int, snip_no::Int, budget::Float64; base_dir::String="data/SNIP/")
    # Define file paths
    intd_arc = joinpath(base_dir, "intd_arc$(instance_no).txt")
    arcgain = joinpath(base_dir, "arcgain$(instance_no).txt")
    scenarios_file = joinpath(base_dir, "Scenarios.txt")
    psi_file = joinpath(base_dir, "psi.txt")

    # Read sensor installation arcs (D)
    D = Vector{Tuple{Int,Int,Float64,Float64}}()
    if isfile(intd_arc)
        for line in eachline(intd_arc)
            line = strip(line)
            isempty(line) && continue
            
            vals = filter(!isempty, split(line, '\t'))

            i, j = parse.(Int, vals[1:2])
            r = parse(Float64, vals[3])
            
            q = if snip_no == 2
                r * 0.5
            elseif snip_no == 3
                r * 0.1
            elseif snip_no == 4
                0.0
            else
                parse(Float64, vals[4])
            end
            
            push!(D, (i, j, r, q))
        end
    end

    # Read non-sensor arcs (A_minus_D)
    A_minus_D = Vector{Tuple{Int,Int,Float64}}()
    if isfile(arcgain)
        for line in eachline(arcgain)
            isempty(strip(line)) && continue
            vals = filter(!isempty, split(line, '\t'))
            i, j = parse.(Int, vals[1:2])
            r = parse(Float64, vals[3])
            push!(A_minus_D, (i, j, r))
        end
    end

    # Read scenarios
    S = Vector{Tuple{Int,Int,Float64}}()
    if isfile(scenarios_file)
        for line in eachline(scenarios_file)
            isempty(strip(line)) && continue
            vals = split(strip(line))
            for i in 1:3:length(vals)
                if i + 2 <= length(vals)
                    try
                        origin = parse(Int, vals[i])
                        dest = parse(Int, vals[i+1])
                        prob = parse(Float64, vals[i+2])
                        push!(S, (origin, dest, prob))
                    catch e
                        @warn "Error processing values at index $i: $(vals[i:i+2])"
                        continue
                    end
                end
            end
        end
    end

    # Read psi matrix
    psi_content = read(psi_file, String)
    psi_content = replace(psi_content, r"\s+" => "")
    arrays = match(r"\[(.*)\]", psi_content).captures[1]
    array_strings = split(arrays, "],[")
    
    psi = Vector{Vector{Float64}}()
    for arr_str in array_strings
        arr_str = replace(arr_str, "[" => "")
        arr_str = replace(arr_str, "]" => "")
        values = split(arr_str, ',')
        row = Float64[]
        for val in values
            val = strip(val)
            if !isempty(val)
                if count(".", val) > 1
                    parts = split(val, '.')
                    val = parts[1] * "." * join(parts[2:end], "")
                end
                push!(row, parse(Float64, val))
            end
        end
        push!(psi, row)
    end

    # Create node mapping and remap indices
    node_mapping = create_node_mapping(D, A_minus_D, S)
    
    # Remap node indices
    D_remapped = [(node_mapping[i], node_mapping[j], r, q) for (i, j, r, q) in D]
    A_minus_D_remapped = [(node_mapping[i], node_mapping[j], r) for (i, j, r) in A_minus_D]
    S_remapped = [(node_mapping[i], node_mapping[j], prob) for (i, j, prob) in S]
    
    return SNIPData(
        length(node_mapping),
        length(S),
        S_remapped,
        D_remapped,
        A_minus_D_remapped,
        psi,
        budget
    )
end

