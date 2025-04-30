using JSON

struct SNIPData <: AbstractData
    num_nodes::Int
    num_scenarios::Int
    scenarios::Vector{Tuple{Int,Int,Float64}} # (from_node, to_node, probability)
    D::Vector{Tuple{Int,Int,Float64,Float64}} # (from_node, to_node, r, q)
    A_minus_D::Vector{Tuple{Int,Int,Float64}} # (from_node, to_node, r)
    ψ::Vector{Vector{Float64}} # ψ matrix
    budget::Float64
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

function read_snip_data(instance_no::Int, snip_no::Int, budget::Float64; base_dir::String="example/snip/data/SNIP/")
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