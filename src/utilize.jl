export assign_solver!

function assign_solver!(model::Model, solver::Symbol)
    if solver == :Gurobi
        set_optimizer(model, Gurobi.Optimizer)
    elseif solver == :CPLEX
        set_optimizer(model, CPLEX.Optimizer)
    else
        error("Unsupported solver: $solver")
    end
    set_optimizer_attribute(model, MOI.Silent(), true)
end

