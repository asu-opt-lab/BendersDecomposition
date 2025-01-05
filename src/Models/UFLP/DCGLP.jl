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

function add_problem_specific_constraints!(model::Model, data::UFLPData, ::AbstractNormType) end