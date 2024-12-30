export DCGLP

mutable struct DCGLP <: AbstractDCGLP
    model::Model
    γ_constraints::Dict{Symbol,Any}
    γ_values::Vector{Tuple{Float64, Vector{Float64}, Union{Float64, Vector{Float64}}}}
    disjunctive_inequalities_constraints::Vector{ConstraintRef}
    dcglp_constraints::Any
    master_cuts::Any
end

# Function to create the DCGLP model
function create_dcglp(data::AbstractData, disjunction_system::DisjunctiveCut)

    model = create_base_model(data, disjunction_system.norm_type)
    add_problem_specific_constraints!(model, data, disjunction_system.norm_type)
    add_t_constraints!(model, data, disjunction_system.base_cut_strategy, disjunction_system.norm_type)
    add_norm_specific_components!(model, data, disjunction_system.base_cut_strategy, disjunction_system.norm_type)

    return DCGLP(model, collect_constraints(model), [], [], [], [])
end

# Function to create the base model for StandardNorm
function create_base_model(data::AbstractData, ::StandardNorm)
    model = Model()
    N = data.n_facilities

    # Define variables
    @variable(model, τ)
    @variable(model, k₀ >= 0)
    @variable(model, kₓ[1:N])
    @variable(model, v₀ >= 0)
    @variable(model, vₓ[1:N])


    # Add constraints
    # @constraint(model, consigma1, τ >= k₀*(constant+1) - coef'kₓ) 
    @constraint(model, coneta1[j in 1:N], τ >= -k₀ + kₓ[j]) 
    # @constraint(model, consigma2, τ >= -v₀*constant + coef'vₓ ) 
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
    N = data.n_facilities

    # Define variables
    @variable(model, τ)
    @variable(model, k₀>=0)
    @variable(model, kₓ[1:N])
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

    # Add γ constraints
    @constraint(model, con0, k₀ + v₀ == 1)
    @constraint(model, conx[i=1:N], kₓ[i] + vₓ[i] - sx[i] == 0)  #x̂

    return model
end

# function add_problem_specific_constraints!(model::Model, data::AbstractData) end

# Generic function to add t constraints
function add_t_constraints!(model::Model, ::AbstractData, ::CutStrategy, ::StandardNorm)
    @variable(model, kₜ)
    @variable(model, vₜ)
    @constraint(model, cont, kₜ + vₜ == 0)
end

function add_t_constraints!(model::Model, ::AbstractData, ::CutStrategy, ::LNorm)
    @variable(model, kₜ)
    @variable(model, vₜ)
    @variable(model, st)
    @constraint(model, cont, kₜ + vₜ - st == 0)
end

# Function to add norm-specific components for StandardNorm
function add_norm_specific_components!(model::Model, data::AbstractData, ::CutStrategy, norm_type::StandardNorm) end

# Function to add norm-specific components for LNorm (ScalarDimension)
function add_norm_specific_components!(model::Model, data::AbstractData, ::CutStrategy, norm_type::LNorm)
    N = data.n_facilities
    if norm_type == L1Norm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.NormInfinityCone(1 + N + 1))
    elseif norm_type == L2Norm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.SecondOrderCone(1 + N + 1))
    elseif norm_type == LInfNorm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.NormOneCone(1 + N + 1))
    else
        error("Unsupported norm type: $(typeof(norm_type))")
    end
end

collect_constraints(model::Model) = Dict{Symbol,Any}(
    :γ₀ => model[:con0],
    :γₓ => model[:conx],
    :γₜ => model[:cont]
)

