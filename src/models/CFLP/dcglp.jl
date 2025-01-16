# ============================================================================
# add_problem_specific_constraints!
# ============================================================================

function add_problem_specific_constraints!(model::Model, data::CFLPData, ::StandardNorm)
    N = data.n_facilities
    @constraint(model, conw1, model[:τ] >= -sum(data.capacities[j]*model[:kₓ][j] for j in 1:N) + sum(data.demands)*model[:k₀])
    @constraint(model, conw2, model[:τ] >= -sum(data.capacities[j]*model[:vₓ][j] for j in 1:N) + sum(data.demands)*model[:v₀])
end

function add_problem_specific_constraints!(model::Model, data::CFLPData, ::LNorm)
    N = data.n_facilities
    @constraint(model, conw1, 0 >= -sum(data.capacities[j]*model[:kₓ][j] for j in 1:N) + sum(data.demands)*model[:k₀])
    @constraint(model, conw2, 0 >= -sum(data.capacities[j]*model[:vₓ][j] for j in 1:N) + sum(data.demands)*model[:v₀])
end

# ============================================================================
# add_t_constraints!
# ============================================================================

function add_t_constraints!(model::Model, ::CFLPData, ::Union{ClassicalCut, KnapsackCut}, ::StandardNorm)
    @variable(model, kₜ)
    @variable(model, vₜ)
    @constraint(model, cont, kₜ + vₜ == 0)
end

function add_t_constraints!(model::Model, ::CFLPData, ::Union{ClassicalCut, KnapsackCut}, ::LNorm)
    @variable(model, kₜ)
    @variable(model, vₜ)
    @variable(model, st)
    @constraint(model, cont, kₜ + vₜ - st == 0)
end

# ============================================================================
# add_norm_specific_components!
# ============================================================================

function add_norm_specific_components!(model::Model, data::CFLPData, ::Union{ClassicalCut, KnapsackCut}, norm_type::LNorm)
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
