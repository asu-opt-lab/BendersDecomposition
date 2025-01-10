export MCNDPMILP, create_milp

struct MCNDPMILP <: AbstractMILP
    model::Model
end

function create_milp(data::MCNDPData)
    model = Model()

    # Set optimizer to silent mode
    set_optimizer_attribute(model, MOI.Silent(), true)

    @variable(model, x[1:data.num_arcs], Bin)  # Binary decision variables for arc selection
    @variable(model, y[1:data.num_arcs, 1:data.num_commodities] >= 0)  # Flow variables
    
    # Objective function
    @objective(model, Min, 
        sum(data.fixed_costs[a] * x[a] for a in 1:data.num_arcs) +
        sum(data.variable_costs[a] * data.demands[k][3] * y[a,k] 
            for a in 1:data.num_arcs, k in 1:data.num_commodities)
    )
    
    # Capacity constraints
    @constraint(model, capacity[a in 1:data.num_arcs],
        sum(data.demands[k][3] * y[a,k] for k in 1:data.num_commodities) <= data.capacities[a] * x[a]
    )
    
    @constraint(model, [a in 1:data.num_arcs, k in 1:data.num_commodities], y[a,k] <= x[a])

    # Flow conservation constraints
    for k in 1:data.num_commodities # for each commodity
        origin, destination, _ = data.demands[k]
        
        for i in 1:data.num_nodes # for each node
            # Calculate inflow and outflow
            inflow = sum(y[a,k] for a in 1:data.num_arcs if data.arcs[a][2] == i)
            outflow = sum(y[a,k] for a in 1:data.num_arcs if data.arcs[a][1] == i)
            
            # Node balance constraint depends on whether it's origin, destination or intermediate
            rhs = if i == origin
                1
            elseif i == destination
                -1
            else
                0
            end
            
            cut = @constraint(model, outflow - inflow == rhs)
        end
    end


    return MCNDPMILP(model)
end
