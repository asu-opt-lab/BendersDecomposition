using BendersDecomposition
using CPLEX, JuMP, Gurobi
using Statistics
# List of r09 instances
i = "10"
instances = ["r$i.$j.dow" for j in 1:9]

# Store solve times
solve_times = Float64[]

for instance in instances
    data = read_mcndp_instance(instance)
    milp = create_milp(data)
    set_optimizer(milp.model, CPLEX.Optimizer)
    
    # Solve and record time
    optimize!(milp.model)
    @info termination_status(milp.model)
    # @info value.(milp.model[:x])
    push!(solve_times, solve_time(milp.model))
    if termination_status(milp.model) == MOI.OPTIMAL
        mip_obj = objective_value(milp.model)
        @info "Instance: $instance, MIP obj: $mip_obj, Solve time: $(solve_time(milp.model))"
    end
end

avg_time = mean(solve_times)
@info "Average solve time across r09 instances: $avg_time seconds"

