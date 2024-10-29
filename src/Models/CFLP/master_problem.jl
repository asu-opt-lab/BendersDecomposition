



mutable struct CFLPMasterProblem <: AbstractMasterProblem
    model::Model
    var::Dict
    coef::Vector{Int}
    x_value::Vector{Float64}
    value_t::Float64
    obj_value::Float64
    data
end

function create_master_problem(::Type{CFLPMasterProblem}, data, solver::Symbol)
    model = if solver == :CPLEX
        Model(CPLEX.Optimizer)
    elseif solver == :Gurobi
        Model(Gurobi.Optimizer)
    else
        error("Unsupported solver: $solver")
    end
    
    set_optimizer_attribute(model, MOI.Silent(), true)

    # Replace println with logging
    @info "Building Master Problem for CFLP"
    
    N = data.n_facilities
    
    @variable(model, x[1:N], Bin)
    @variable(model, t >= -1e6)

    @objective(model, Min, sum(data.fixed_costs[i] * x[i] for i in 1:N) + t)

    # Ensure total capacity meets total demand
    @constraint(model, sum(data.capacities[i] * x[i] for i in 1:N) >= sum(data.demands))

    # Initialize the master problem struct
    mp = CFLPMasterProblem(
        model,
        Dict("x" => x, "t" => t),  
        data.fixed_costs,
        zeros(N),
        0.0,
        0.0,
        data
    )

    return mp
end


