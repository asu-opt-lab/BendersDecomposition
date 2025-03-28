# ============================================================================
# add_problem_specific_constraints!
# ============================================================================

function add_problem_specific_constraints!(model::Model, data::UFLPData, ::LNorm)
    @constraint(model, 0 >= - sum(model[:kₓ]) + 2*model[:k₀])
    @constraint(model, 0 >= - sum(model[:vₓ]) + 2*model[:v₀])
end
function add_problem_specific_constraints!(model::Model, data::UFLPData, ::StandardNorm)
    @constraint(model, model[:τ] >= - sum(model[:kₓ]) + 2*model[:k₀])
    @constraint(model, model[:τ] >= - sum(model[:vₓ]) + 2*model[:v₀])
end

# ============================================================================
# add_t_constraints!
# ============================================================================
function add_t_constraints!(model::Model, ::UFLPData, ::ClassicalCut, ::StandardNorm)
    @variable(model, kₜ)
    @variable(model, vₜ)
    @constraint(model, cont, kₜ + vₜ == 0)
end

function add_t_constraints!(model::Model, ::UFLPData, ::ClassicalCut, ::LNorm)
    @variable(model, kₜ)
    @variable(model, vₜ)
    @constraint(model, kₜ >= -1e6 * model[:k₀])
    @constraint(model, vₜ >= -1e6 * model[:v₀])
    @variable(model, st)
    @constraint(model, cont, kₜ + vₜ - st == 0)
end

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
    @constraint(model, kₜ .>= -1e6 * model[:k₀])
    @constraint(model, vₜ .>= -1e6 * model[:v₀])
    @variable(model, st[1:M])
    γₜconstarint = @constraint(model, cont[i=1:M], kₜ[i] + vₜ[i] - st[i] == 0)
    return γₜconstarint
end

# ============================================================================
# add_norm_specific_components!
# ============================================================================

function add_norm_specific_components!(model::Model, data::UFLPData, ::FatKnapsackCut, norm_type::LNorm)
    N = data.n_facilities
    M = data.n_customers
    dim = 1 + N + M
    if norm_type == L1Norm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.NormInfinityCone(dim))
    elseif norm_type == L2Norm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.SecondOrderCone(dim))
    elseif norm_type == LInfNorm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.NormOneCone(dim))
    else
        error("Unsupported norm type: $(typeof(norm_type))")
    end
end

function add_norm_specific_components!(model::Model, data::UFLPData, ::ClassicalCut, norm_type::LNorm)
    N = data.n_facilities
    dim = 1 + N + 1
    if norm_type == L1Norm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.NormInfinityCone(dim))
    elseif norm_type == L2Norm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.SecondOrderCone(dim))
    elseif norm_type == LInfNorm()
        @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.NormOneCone(dim))
    else
        error("Unsupported norm type: $(typeof(norm_type))")
    end
end