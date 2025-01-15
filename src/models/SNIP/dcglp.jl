# ============================================================================
# add_problem_specific_constraints!
# ============================================================================

function add_problem_specific_constraints!(model::Model, data::SNIPData, ::AbstractNormType)
end

# ============================================================================
# add_t_constraints!
# ============================================================================

function add_t_constraints!(model::Model, data::SNIPData, ::ClassicalCut, ::StandardNorm)
    M = data.num_scenarios
    @variable(model, kₜ[1:M])
    @variable(model, vₜ[1:M])
    γₜconstarint = @constraint(model, cont[i=1:M], kₜ[i] + vₜ[i] == 0)
    return γₜconstarint
end

function add_t_constraints!(model::Model, data::SNIPData, ::ClassicalCut, ::LNorm)
    
    M = data.num_scenarios
    @variable(model, kₜ[1:M])
    @variable(model, vₜ[1:M])
    @variable(model, st[1:M])
    γₜconstarint = @constraint(model, cont[i=1:M], kₜ[i] + vₜ[i] - st[i] == 0)
    return γₜconstarint
end

# ============================================================================
# add_norm_specific_components!
# ============================================================================

function add_norm_specific_components!(model::Model, data::SNIPData, ::ClassicalCut, norm_type::LNorm)
    N = length(data.D)
    M = data.num_scenarios
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