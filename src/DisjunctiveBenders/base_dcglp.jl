export DCGLP

# Function to create the DCGLP model
function create_dcglp(master::AbstractMasterProblem, disjunction_system::DisjunctiveCut)

    model = create_base_model(master, disjunction_system.norm_type)
    # add_problem_specific_constraints!(model, master, disjunction_system.norm_type)
    add_t_constraints!(model, master, disjunction_system.norm_type)
    add_norm_specific_components!(model, master, disjunction_system.norm_type)

    return DCGLP(model, collect_constraints(model), [], [], ConstraintRef[], AffExpr[])
end

# ============================================================================
# create_base_model
# ============================================================================
function create_base_model(master::AbstractMasterProblem, ::StandardNorm)
    model = Model()

    N = length(master.integer_variable_values)

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
function create_base_model(master::AbstractMasterProblem, ::LNorm)
    # Create model with specified solver
    model = Model()

    # Define problem dimensions
    N = length(master.integer_variable_values)

    # Define variables
    @variable(model, τ)
    @variable(model, k₀)
    @variable(model, kₓ[1:N])
    @variable(model, v₀)
    @variable(model, vₓ[1:N])
    @variable(model, sₓ[1:N])

    # Set objective
    @objective(model, Min, τ)

    # Add constraints
    @constraint(model, coneta1[j in 1:N], 0 >= -k₀ + kₓ[j]) 
    @constraint(model, coneta2[j in 1:N], 0 >= -v₀ + vₓ[j])
    @constraint(model, conv1[j in 1:N], 0 >= -kₓ[j])
    @constraint(model, conv2[j in 1:N], 0 >= -vₓ[j])
    @constraint(model, conk0, k₀ >= 0)
    @constraint(model, conv0, v₀ >= 0)

    # Add γ constraints
    @constraint(model, con0, k₀ + v₀ == 1)
    @constraint(model, conx[i=1:N], kₓ[i] + vₓ[i] - sₓ[i] == 0)  #x̂

    return model
end

# ============================================================================
# add_problem_specific_constraints!
# ============================================================================
function add_problem_specific_constraints!(model::Model, master::AbstractMasterProblem) 
    @warn "No problem specific constraints implemented"
end

# ============================================================================
# add_t_constraints!
# ============================================================================

function add_t_constraints!(model::Model, master::GenericMasterProblem, norm_type::LNorm)
    @variable(model, kₜ)
    @variable(model, vₜ)
    @variable(model, sₜ)
    @constraint(model, cont, kₜ + vₜ - sₜ == 0)
end

# ============================================================================
# add_norm_specific_components!
# ============================================================================

function add_norm_specific_components!(model::Model, master::GenericMasterProblem, norm_type::LNorm)
    N = length(master.integer_variable_values)
    dim = 1 + N + 1
    if norm_type == L1Norm()
        @constraint(model, concone, [model[:τ]; model[:sₓ]; model[:sₜ]] in MOI.NormInfinityCone(dim))
    elseif norm_type == L2Norm()
        @constraint(model, concone, [model[:τ]; model[:sₓ]; model[:sₜ]] in MOI.SecondOrderCone(dim))
    elseif norm_type == LInfNorm()
        @constraint(model, concone, [model[:τ]; model[:sₓ]; model[:sₜ]] in MOI.NormOneCone(dim))
    else
        error("Unsupported norm type: $(typeof(norm_type))")
    end
end

# ============================================================================
# helper functions
# ============================================================================
collect_constraints(model::Model) = Dict{Symbol,Any}(
    :γ₀ => model[:con0],
    :γₓ => model[:conx],
    :γₜ => model[:cont]
)

