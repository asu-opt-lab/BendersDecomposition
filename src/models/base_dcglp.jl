export DCGLP

"""
    DCGLP <: AbstractDCGLP

A mutable struct representing a Disjunctive Cut Generation Linear Program (DCGLP).

# Fields
- `model::Model`: The underlying JuMP optimization model
- `γ_constraints::Dict{Symbol,Any}`: Dictionary storing the gamma constraints of the model (γ₀, γₓ, γₜ)
- `γ_values::Vector{Tuple{Float64, Vector{Float64}, Union{Float64, Vector{Float64}}}}`: Vector storing gamma values
- `disjunctive_inequalities_constraints::Vector{ConstraintRef}`: Vector of disjunctive inequality constraints for each index selected
- `dcglp_constraints::Any`: Storage for disjunctive cuts
- `master_cuts::Any`: Storage for cuts being added to the master problem

The DCGLP is used to generate cutting planes for disjunctive programming problems,
particularly useful in mixed-integer programming.
"""
mutable struct DCGLP <: AbstractDCGLP
    model::Model
    var_kv_x::Dict
    γ_constraints::Dict{Symbol,Any}
    γ_values::Vector{Tuple{Float64, Vector{Float64}, Union{Float64, Vector{Float64}}}}
    disjunctive_inequalities_constraints::Vector{ConstraintRef}
    dcglp_constraints::Any
    master_cuts::Any
    strengthen_used::Vector{Int64}
    disjuncindex_cut_pairs::Dict
    prev_dcglp_disjuncuts::Any
    split_index::Vector{Int64}
end

# Function to create the DCGLP model
function create_dcglp(data::AbstractData, disjunction_system::DisjunctiveCut)

    model, var_kv_x = create_base_model(data, disjunction_system.norm_type)
    add_problem_specific_constraints!(model, data, disjunction_system.norm_type)
    add_t_constraints!(model, data, disjunction_system.base_cut_strategy, disjunction_system.norm_type)
    add_norm_specific_components!(model, data, disjunction_system.base_cut_strategy, disjunction_system.norm_type)

    return DCGLP(model, var_kv_x, collect_constraints(model), [], [], [], [], [], Dict(), [], [])
end

# ============================================================================
# create_base_model
# ============================================================================
function create_base_model(data::AbstractData, ::StandardNorm)
    model = Model()

    N = get_problem_size(data)

    # Define variables
    @variable(model, τ)
    @variable(model, k₀ >= 0)
    @variable(model, kₓ[1:N])
    @variable(model, v₀ >= 0)
    @variable(model, vₓ[1:N])


    # Add constraints
    @constraint(model, coneta1[j in 1:N], τ >= -k₀ + kₓ[j]) 
    @constraint(model, coneta2[j in 1:N], τ >= -v₀ + vₓ[j])
    @constraint(model, conv1[j in 1:N], τ >= -kₓ[j])
    @constraint(model, conv2[j in 1:N], τ >= -vₓ[j])

    # Add γ constraints
    @constraint(model, con0, k₀ + v₀ == 1)
    @constraint(model, conx[i=1:N], kₓ[i] + vₓ[i] == 0)  #-x̂

    return model
end

# Function to create the base model for LNorm
function create_base_model(data::AbstractData, ::LNorm)
    # Create model with specified solver
    model = Model()

    # Define problem dimensions
    N = get_problem_size(data)

    # Define variables
    @variable(model, τ)
    # @variable(model, k₀)
    @variable(model, k₀>=0)
    @variable(model, kₓ[1:N])
    # @variable(model, v₀)
    @variable(model, v₀>=0)
    @variable(model, vₓ[1:N])
    @variable(model, sx[1:N])

    # Set objective
    @objective(model, Min, τ)

    # Add constraints
    @constraint(model, coneta1[j in 1:N], 0 >= -k₀ + kₓ[j]) 
    @constraint(model, coneta2[j in 1:N], 0 >= -v₀ + vₓ[j])
    @constraint(model, conv1[j in 1:N], 0 >= -kₓ[j])
    @constraint(model, conv2[j in 1:N], 0 >= -vₓ[j])
    # @constraint(model, conk0, k₀>=0)
    # @constraint(model, conv0, v₀>=0)

    # Add γ constraints
    @constraint(model, con0, k₀ + v₀ == 1)
    @constraint(model, conx[i=1:N], kₓ[i] + vₓ[i] - sx[i] == 0)  #x̂

    return model, Dict(:kₓ => kₓ, :vₓ => vₓ)
end

# ============================================================================
# add_problem_specific_constraints!
# ============================================================================
function add_problem_specific_constraints!(model::Model, data::AbstractData) 
    @warn "No problem specific constraints implemented for data type: $(typeof(data))"
end

# ============================================================================
# add_t_constraints!
# ============================================================================

function add_t_constraints!(model::Model, data::AbstractData, cut_strategy::CutStrategy, norm_type::AbstractNormType)
    @error "No t constraints implemented for data type: $(typeof(data)) and cut strategy: $(typeof(cut_strategy)) and norm type: $(typeof(norm_type))"
end

# ============================================================================
# add_norm_specific_components!
# ============================================================================

function add_norm_specific_components!(model::Model, data::AbstractData, cut_strategy::CutStrategy, norm_type::StandardNorm) end


function add_norm_specific_components!(model::Model, data::AbstractData, cut_strategy::CutStrategy, norm_type::LNorm)
    @error "No norm specific components implemented for data type: $(typeof(data)) and cut strategy: $(typeof(cut_strategy)) and norm type: $(typeof(norm_type))"
end

# ============================================================================
# helper functions
# ============================================================================
collect_constraints(model::Model) = Dict{Symbol,Any}(
    :γ₀ => model[:con0],
    :γₓ => model[:conx],
    :γₜ => model[:cont]
)

get_problem_size(data::SNIPData) = length(data.D)
get_problem_size(data::Union{CFLPData, UFLPData, SCFLPData, MCNDPData}) = data.n_facilities
get_problem_size(data::Any) = error("Unsupported data type: $(typeof(data)). " *
    "Supported types are: SNIPData, CFLPData, UFLPData, SCFLPData, MCNDPData")