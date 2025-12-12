export UFLKnapsackOracle, generate_cuts, calculate_KP_value
mutable struct UFLKnapsackOracleParam <: AbstractOracleParam
    slim::Bool
    add_only_violated_cuts::Bool
    rtol::Float64

    function UFLKnapsackOracleParam(; slim = false, add_only_violated_cuts = false, rtol = 1e-9)
        new(slim, add_only_violated_cuts, rtol)
    end
end
mutable struct UFLKnapsackOracle <: AbstractTypicalOracle
    
    param::UFLKnapsackOracleParam
    
    sorted_cost_demands::Vector{Vector{Float64}}
    sorted_indices::Vector{Vector{Int}}

    J::Int
    obj_values::Vector{Float64}

    function UFLKnapsackOracle(data::UFLPData; 
        scen_idx::Int=-1, 
        param::UFLKnapsackOracleParam = UFLKnapsackOracleParam())
            @debug "Building knapsack oracle for UFLP"
            
            J = data.n_customers
            cost_demands = [data.costs[:,j] .* data.demands[j] for j in 1:J]
            sorted_indices = [sortperm(cost_demands[j]) for j in 1:J]
            sorted_cost_demands = [cost_demands[j][sorted_indices[j]] for j in 1:J]

            obj_values = Vector{Float64}(undef, J)

            new(param, sorted_cost_demands, sorted_indices, J, obj_values)
    end

    UFLKnapsackOracle() = new()
end

function generate_cuts(oracle::UFLKnapsackOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600.0)
    tic = time()
    critical_facility = Vector{Int}(undef, oracle.J)
    for j in 1:oracle.J
        sorted_indices = oracle.sorted_indices[j]
        c_sorted = oracle.sorted_cost_demands[j]
        x_sorted = x_value[sorted_indices]

        # Find critical item and calculate contribution
        k = find_critical_item(c_sorted, x_sorted)

        # Calculate objective value contribution
        oracle.obj_values[j] = c_sorted[k] - (k > 1 ? sum((c_sorted[k] - c_sorted[i]) * x_sorted[i] for i in 1:k-1) : 0)

        if oracle.obj_values[j] >= t_value[j] * (1 + oracle.param.rtol)
            critical_facility[j] = k
        else
            critical_facility[j] = oracle.param.add_only_violated_cuts ? -1 : k
        end
    end
    
    if get_sec_remaining(tic, time_limit) <= 0.0
        throw(TimeLimitException("Time limit reached during cut generation"))
    end

    customers = findall(x -> x != -1, critical_facility)

    # is_in_L should be determined by the sum of t's, must not individually
    is_in_L = sum(oracle.obj_values) >= sum(t_value) * (1 + oracle.param.rtol) ? false : true

    hyperplanes = Vector{Hyperplane}()
    for j in customers
        k = critical_facility[j] 
        sorted_indices = oracle.sorted_indices[j]
        c_sorted = oracle.sorted_cost_demands[j]
        
        h = Hyperplane(length(x_value), oracle.J)
        h.a_t[j] = -1.0
        h.a_0 = c_sorted[k]
        for i=1:k-1
            h.a_x[sorted_indices[i]] = -(c_sorted[k] - c_sorted[i])
        end
        push!(hyperplanes, h)
    end

    if is_in_L
        return !(oracle.param.slim) ? (true, hyperplanes, deepcopy(t_value)) : (true, [aggregate(hyperplanes)], deepcopy(t_value))
    end
    return !(oracle.param.slim) ? (false, hyperplanes, deepcopy(oracle.obj_values)) : (false, [aggregate(hyperplanes)], deepcopy(oracle.obj_values))
end

function find_critical_item(c::Vector{Float64}, x::Vector{Float64})
    
    sum_x::Float64 = 0.0
    for (idx, val) in enumerate(x)
        sum_x += val
        if sum_x >= 1.0
            return idx
        end
    end
    throw(AlgorithmException("`k` cannot be `nothing` as sum(x) >= 2 is enforced. Check the models."))
end