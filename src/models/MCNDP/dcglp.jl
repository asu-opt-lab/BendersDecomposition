##### NOT IMPLEMENTED YET #####

# ============================================================================
# add_problem_specific_constraints!
# ============================================================================

function add_problem_specific_constraints!(model::Model, data::MCNDPData, ::AbstractNormType) end

# ============================================================================
# add_t_constraints!
# ============================================================================

# function add_t_constraints!(model::Model, ::CFLPData, ::Union{ClassicalCut, KnapsackCut}, ::StandardNorm)
#     @variable(model, kₜ)
#     @variable(model, vₜ)
#     @constraint(model, cont, kₜ + vₜ == 0)
# end

# function add_t_constraints!(model::Model, ::CFLPData, ::Union{ClassicalCut, KnapsackCut}, ::LNorm)
#     @variable(model, kₜ)
#     @variable(model, vₜ)
#     @variable(model, st)
#     @constraint(model, cont, kₜ + vₜ - st == 0)
# end

# ============================================================================
# add_norm_specific_components!
# ============================================================================

# function add_norm_specific_components!(model::Model, data::CFLPData, ::Union{ClassicalCut, KnapsackCut}, norm_type::LNorm)
#     N = data.n_facilities
#     if norm_type == L1Norm()
#         @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.NormInfinityCone(1 + N + 1))
#     elseif norm_type == L2Norm()
#         @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.SecondOrderCone(1 + N + 1))
#     elseif norm_type == LInfNorm()
#         @constraint(model, concone, [model[:τ]; model[:sx]; model[:st]] in MOI.NormOneCone(1 + N + 1))
#     else
#         error("Unsupported norm type: $(typeof(norm_type))")
#     end
# end

