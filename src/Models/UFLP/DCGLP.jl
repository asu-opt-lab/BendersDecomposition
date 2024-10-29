export UFLPDCGLP
# Mutable struct representing the DCGLP (Disjunctive Cut Generating Linear Program) for UFLP
mutable struct UFLPDCGLP <: AbstractDCGLP
    model::Model
    γconstraints::Dict{Symbol,Any}
end

# Function to create the DCGLP model
function create_dcglp(
    data::UFLPData,
    disjunctive_inequality::Tuple{Vector{Int}, Int},
    _cut_strategy::CutGenerationStrategy,
    norm_type::AbstractNormType
)
    # Create base model
    model,γ₀constarint,γₓconstarint = create_base_model(data, disjunctive_inequality,norm_type)
    # Add t constraints
    γₜconstarint = add_t_constraints!(model, data, _cut_strategy, norm_type)
    # Add norm-specific components
    add_norm_specific_components!(model, data, _cut_strategy, norm_type)
    # Return UFLPDCGLP struct
    return UFLPDCGLP(model, Dict(:γ₀=>γ₀constarint,:γₓ=>γₓconstarint,:γₜ=>γₜconstarint))
end

# Function to create the base model for StandardNorm
function create_base_model(data::UFLPData, disjunctive_inequality::Tuple{Vector{Int}, Int}, ::StandardNorm)

    coef,constant = disjunctive_inequality

    model = Model()

    # Set optimizer to silent mode
    set_optimizer_attribute(model, MOI.Silent(), true)

    # Define problem dimensions
    N = data.n_facilities

    # Define variables
    @variable(model, τ)
    @variable(model, k₀ >= 0)
    @variable(model, kₓ[1:N])
    @variable(model, v₀ >= 0)
    @variable(model, vₓ[1:N])

    total_demands = sum(data.demands)

    # Add constraints
    @constraint(model, consigma1, τ >= k₀*(constant+1) - coef'kₓ) 
    @constraint(model, coneta1[j in 1:N], τ >= -k₀ + kₓ[j]) 
    @constraint(model, consigma2, τ >= -v₀*constant + coef'vₓ ) 
    @constraint(model, coneta2[j in 1:N], τ >= -v₀ + vₓ[j])
    @constraint(model, conv1[j in 1:N], τ >= -kₓ[j])
    @constraint(model, conv2[j in 1:N], τ >= -vₓ[j])

    # Add γ constraints
    γ₀constarint = @constraint(model, con0, k₀ + v₀ == 1)
    γₓconstarint = @constraint(model, conx[i=1:N], kₓ[i] + vₓ[i] == 0)  #-x̂

    return model,γ₀constarint,γₓconstarint
end

# Function to create the base model for LNorm
function create_base_model(data::UFLPData, disjunctive_inequality::Tuple{Vector{Int}, Int}, ::LNorm)

    coef,constant = disjunctive_inequality

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
    @variable(model, sx[i = 1:N])

    # Set objective
    @objective(model, Min, τ)

    # Add constraints
    @constraint(model, consigma1, 0 >= k₀*(constant+1) - coef'kₓ) 
    @constraint(model, coneta1[j in 1:N], 0 >= -k₀ + kₓ[j]) 
    @constraint(model, consigma2, 0 >= -v₀*constant + coef'vₓ) 
    @constraint(model, coneta2[j in 1:N], 0 >= -v₀ + vₓ[j])
    @constraint(model, conv1[j in 1:N], 0 >= -kₓ[j])
    @constraint(model, conv2[j in 1:N], 0 >= -vₓ[j])

    # Add γ constraints
    γ₀constarint = @constraint(model, con0, k₀ + v₀ == 1)
    γₓconstarint = @constraint(model, conx[i=1:N], kₓ[i] + vₓ[i] - sx[i] == 0)  #x̂

    return model,γ₀constarint,γₓconstarint
end

# Generic function to add t constraints
function add_t_constraints!(model::Model, ::UFLPData, ::CutGenerationStrategy, ::StandardNorm)
    
    @variable(model, kₜ)
    @variable(model, vₜ)
    γₜconstarint = @constraint(model, cont, kₜ + vₜ == 0)
    return γₜconstarint
end

function add_t_constraints!(model::Model, ::UFLPData, ::CutGenerationStrategy, ::LNorm)
    @variable(model, kₜ)
    @variable(model, vₜ)
    @variable(model, st)
    γₜconstarint = @constraint(model, cont, kₜ + vₜ - st == 0)
    return γₜconstarint
end

# for multiple t variables
function add_t_constraints!(model::Model, data::UFLPData, ::FatKnapsackCut, ::StandardNorm)
    M = data.n_customers
    @variable(model, kₜ[1:M])
    @variable(model, vₜ[1:M])
    γₜconstarint = @constraint(model, cont[i=1:M], kₜ[i] + vₜ[i] == 0)
    return γₜconstarint
end

function add_t_constraints!(model::Model, data::UFLPData, ::FatKnapsackCut, ::LNorm)
    
    M = data.n_customers
    @variable(model, kₜ[1:M])
    @variable(model, vₜ[1:M])
    @variable(model, st[1:M])
    γₜconstarint = @constraint(model, cont[i=1:M], kₜ[i] + vₜ[i] - st[i] == 0)
    return γₜconstarint
end

# Function to add norm-specific components for StandardNorm
function add_norm_specific_components!(model::Model, data::UFLPData, ::CutGenerationStrategy, norm_type::StandardNorm)
end

# Function to add norm-specific components for LNorm (ScalarDimension)
function add_norm_specific_components!(model::Model, data::UFLPData, ::CutGenerationStrategy, norm_type::LNorm)
    N = data.n_facilities
    if norm_type == L1Norm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.NormInfinityCone(1 + N + 1))
    elseif norm_type == L2Norm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.SecondOrderCone(1 + N + 1))
    elseif norm_type == LInfNorm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.NormInfinityCone(1 + N + 1))
    else
        error("Unsupported norm type: $(typeof(norm_type))")
    end
end

# Function to add norm-specific components for FatKnapsackCut
function add_norm_specific_components!(model::Model, data::UFLPData, ::FatKnapsackCut, norm_type::LNorm)
    N = data.n_facilities
    M = data.n_customers
    if norm_type == L1Norm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.NormInfinityCone(1 + N + M))
    elseif norm_type == L2Norm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.SecondOrderCone(1 + N + M))
    elseif norm_type == LInfNorm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.NormInfinityCone(1 + N + M))
    else
        error("Unsupported norm type: $(typeof(norm_type))")
    end
end

