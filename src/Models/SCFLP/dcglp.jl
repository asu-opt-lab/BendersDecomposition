function add_problem_specific_constraints!(model::Model, data::SCFLPData, ::StandardNorm)
    N = data.n_facilities
    S = data.n_scenarios
    @constraint(model, conw1[s = 1:S], model[:τ] >= -sum(data.capacities[j]*model[:kₓ][j] for j in 1:N) + sum(data.demands[s])*model[:k₀])
    @constraint(model, conw2[s = 1:S], model[:τ] >= -sum(data.capacities[j]*model[:vₓ][j] for j in 1:N) + sum(data.demands[s])*model[:v₀])
end

function add_problem_specific_constraints!(model::Model, data::SCFLPData, ::LNorm)
    N = data.n_facilities
    S = data.n_scenarios
    @constraint(model, conw1[s = 1:S], 0 >= -sum(data.capacities[j]*model[:kₓ][j] for j in 1:N) + sum(data.demands[s])*model[:k₀])
    @constraint(model, conw2[s = 1:S], 0 >= -sum(data.capacities[j]*model[:vₓ][j] for j in 1:N) + sum(data.demands[s])*model[:v₀])
end

# for multiple t variables
function add_t_constraints!(model::Model, data::SCFLPData, ::Union{ClassicalCut, KnapsackCut}, ::StandardNorm)
    M = data.n_scenarios
    @variable(model, kₜ[1:M])
    @variable(model, vₜ[1:M])
    γₜconstarint = @constraint(model, cont[i=1:M], kₜ[i] + vₜ[i] == 0)
    return γₜconstarint
end

function add_t_constraints!(model::Model, data::SCFLPData, ::Union{ClassicalCut, KnapsackCut}, ::LNorm)
    
    M = data.n_scenarios
    @variable(model, kₜ[1:M])
    @variable(model, vₜ[1:M])
    @variable(model, st[1:M])
    γₜconstarint = @constraint(model, cont[i=1:M], kₜ[i] + vₜ[i] - st[i] == 0)
    return γₜconstarint
end


function add_norm_specific_components!(model::Model, data::SCFLPData, ::Union{ClassicalCut, KnapsackCut}, norm_type::LNorm)
    N = data.n_facilities
    M = data.n_scenarios
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