using JSON

struct MCNDPData <: AbstractData
    num_nodes::Int      # Number of nodes
    num_arcs::Int       # Number of arcs
    num_commodities::Int # Number of commodities
    arcs::Vector{Tuple{Int,Int}}  # Arcs (from_node, to_node)
    fixed_costs::Vector{Float64}   # Fixed costs for opening arcs
    variable_costs::Vector{Float64} # Variable costs per unit flow
    capacities::Vector{Float64}     # Arc capacities
    demands::Vector{Tuple{Int,Int,Float64}} # Demands (origin, destination, quantity)
end

function read_mcndp_instance(filename::String;filepath="example/mcndp/data/NDR/"::AbstractString)
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