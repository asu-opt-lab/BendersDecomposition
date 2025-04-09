mutable struct UFLKnapsackOracle <: AbstractTypicalOracle
    sorted_cost_demands::Vector{Vector{Float64}}
    sorted_indices::Vector{Vector{Int}}

    J::Int
    obj_values::Vector{Float64}

    slim::Bool
    add_only_violated_cuts::Bool

    function UFLKnapsackOracle(data; slim=false, scen_idx::Int=-1, add_only_violated_cuts=false)
        @debug "Building knapsack oracle for UFLP"
        J = data.problem.n_customers
        cost_demands = [data.problem.costs[:,j] .* data.problem.demands[j] for j in 1:J]
        sorted_indices = [sortperm(cost_demands[j]) for j in 1:J]
        sorted_cost_demands = [cost_demands[j][sorted_indices[j]] for j in 1:J]

        obj_values = Vector{Float64}(undef, J)

        new(sorted_cost_demands, sorted_indices, J, obj_values, slim, add_only_violated_cuts)
    end

    UFLKnapsackOracle() = new()
end

function generate_cuts(oracle::UFLKnapsackOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol = 1e-6)
    critical_facility = Vector{Int}(undef, oracle.J)
    for j in 1:oracle.J
        sorted_indices = oracle.sorted_indices[j]
        c_sorted = oracle.sorted_cost_demands[j]
        x_sorted = x_value[sorted_indices]

        # Find critical item and calculate contribution
        k = find_critical_item(c_sorted, x_sorted)

        # Calculate objective value contribution
        oracle.obj_values[j] = c_sorted[k] - (k > 1 ? sum((c_sorted[k] - c_sorted[i]) * x_sorted[i] for i in 1:k-1) : 0)

        if oracle.obj_values[j] >= t_value[j] + tol
            critical_facility[j] = k
        else
            critical_facility[j] = oracle.add_only_violated_cuts ? -1 : k
        end
    end

    # is_in_L should be determined by the sum of t's, must not individually
    is_in_L = sum(oracle.obj_values) >= sum(t_value) + tol ? false : true

    customers = findall(x -> x != -1, critical_facility)
    
    if is_in_L
        return true, [Hyperplane(length(x_value), oracle.J)], t_value
    end

    hyperplanes = Vector{Hyperplane}()
    
    for j in customers
        k = critical_facility[j] 
        sorted_indices = oracle.sorted_indices[j]
        c_sorted = oracle.sorted_cost_demands[j]

        h = Hyperplane(length(x_value), oracle.J)
        h.a_t[j] = -1
        h.a_0 = c_sorted[k]
        for i=1:k-1
            h.a_x[sorted_indices[i]] = -(c_sorted[k] - c_sorted[i])
        end
        push!(hyperplanes, h)
    end
    return !(oracle.slim) ? (false, hyperplanes, oracle.obj_values) : (false, [aggregate(hyperplanes)], oracle.obj_values)
end

function find_critical_item(c::Vector{Float64}, x::Vector{Float64})
    cumsum_x = cumsum(x)
    k = findfirst(>=(1.0), cumsum_x)
    return k === nothing ? length(c) : k
end