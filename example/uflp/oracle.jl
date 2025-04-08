mutable struct UFLKnapsackOracle <: AbstractTypicalOracle
    sorted_cost_demands::Vector{Vector{Float64}}
    sorted_indices::Vector{Vector{Int}}
    selected_k::Dict

    J::Int
    critical_pairs::Vector{Tuple{Int,Int}} # (index of the facility, critical item)
    obj_values::Vector{Float64}

    slim::Bool

    function UFLKnapsackOracle(data; slim=false, scen_idx::Int=-1)
        @debug "Building knapsack oracle for UFLP"
        J = data.problem.n_customers
        cost_demands = [data.problem.costs[:,j] .* data.problem.demands[j] for j in 1:J]
        sorted_indices = [sortperm(cost_demands[j]) for j in 1:J]
        sorted_cost_demands = [cost_demands[j][sorted_indices[j]] for j in 1:J]
        selected_k = Dict(j => [] for j in 1:J)

        critical_pairs = Vector{Tuple{Int,Int}}(undef, J) # (index of the facility, critical item)
        obj_values = Vector{Float64}(undef, J)

        new(sorted_cost_demands, sorted_indices, selected_k, J, critical_pairs, obj_values, slim)
    end

    UFLKnapsackOracle() = new()
end

function generate_cuts(oracle::UFLKnapsackOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol = 1e-6)
    # Process each facility
    # is_violated = false
    for j in 1:oracle.J
        sorted_indices = oracle.sorted_indices[j]
        c_sorted = oracle.sorted_cost_demands[j]
        x_sorted = x_value[sorted_indices]

        # Find critical item and calculate contribution
        k = find_critical_item(c_sorted, x_sorted)

        # Calculate objective value contribution
        oracle.obj_values[j] = c_sorted[k] - (k > 1 ? sum((c_sorted[k] - c_sorted[i]) * x_sorted[i] for i in 1:k-1) : 0)

        if oracle.obj_values[j] >= t_value[j] + tol
            oracle.critical_pairs[j] = (j, k)
            # is_violated = true
            # @info oracle.obj_values[j], t_value[j] 
        else
            oracle.critical_pairs[j] = (j, -1)
            # oracle.critical_pairs[j] = (j, k)
        end
    end

    violated = filter(x -> x[2] != -1, oracle.critical_pairs)
    if length(violated) == 0
        return true, ([zeros(length(t_value))], [zeros(length(x_value))], [0.0]), t_value
    end
    # if !is_violated #length(violated) == 0
    #     return true, ([zeros(length(t_value))], [zeros(length(x_value))], [0.0]), t_value
    # end

    hyperplanes = Vector{Hyperplane}()
    
    # for (index, critical_item) in oracle.critical_pairs
    for (index, critical_item) in violated
        k = critical_item
        sorted_indices = oracle.sorted_indices[index]
        c_sorted = oracle.sorted_cost_demands[index]

        h = Hyperplane(length(x_value), oracle.J)
        h.a_t[index] = -1
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
    k = findfirst(>=(1), cumsum_x)
    return k === nothing ? length(c) : k
end