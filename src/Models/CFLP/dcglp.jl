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

