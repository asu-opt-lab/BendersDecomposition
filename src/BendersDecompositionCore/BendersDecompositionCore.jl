export create_master_problem, create_sub_problem

export CFLPData, UFLPData, SCFLPData, MCNDPData,SNIPData
abstract type AbstractData end

struct CFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    capacities::Vector{Float64}
    demands::Vector{Float64}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

struct UFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    demands::Vector{Float64}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

struct SCFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    num_scenarios::Int
    capacities::Vector{Float64}
    demands::Vector{Vector{Float64}}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

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

struct SNIPData <: AbstractData
    num_nodes::Int
    num_scenarios::Int
    scenarios::Vector{Tuple{Int,Int,Float64}} # (from_node, to_node, probability)
    D::Vector{Tuple{Int,Int,Float64,Float64}} # (from_node, to_node, r, q)
    A_minus_D::Vector{Tuple{Int,Int,Float64}} # (from_node, to_node, r)
    ψ::Vector{Vector{Float64}} # ψ matrix
    budget::Float64
end


struct KnapsackCut <: AbstractCutStrategy end
struct FatKnapsackCut <: AbstractCutStrategy end
struct SlimKnapsackCut <: AbstractCutStrategy end


function create_master_problem end
function create_sub_problem end
function create_master_problem(data::AbstractData, cut_strategy::AbstractCutStrategy)
    throw(MethodError(create_master_problem, (data, cut_strategy)))
end
function create_sub_problem(data::AbstractData, cut_strategy::AbstractCutStrategy)
    throw(MethodError(create_sub_problem, (data, cut_strategy)))
end

include("models/models.jl")
include("utils/utils.jl")