mutable struct UFLKnapsackOracleParam <: AbstractOracleParam
    slim::Bool
    add_only_violated_cuts::Bool

    function UFLKnapsackOracleParam(; slim = false, add_only_violated_cuts = false)
        new(slim, add_only_violated_cuts)
    end
end
mutable struct UFLKnapsackOracle <: AbstractTypicalOracle
    
    oracle_param::UFLKnapsackOracleParam
    
    sorted_cost_demands::Vector{Vector{Float64}}
    sorted_indices::Vector{Vector{Int}}

    J::Int
    obj_values::Vector{Float64}

    function UFLKnapsackOracle(data; 
                               scen_idx::Int=-1, 
                               solver_param::Any = nothing,
                               oracle_param::UFLKnapsackOracleParam = UFLKnapsackOracleParam())
        @debug "Building knapsack oracle for UFLP"
        if solver_param != nothing
            throw(AlgorithmException("UFLKnapsackOracle is solver-free, but it is given with solver_param of type $(type(solver_param))"))
        end
        
        J = data.problem.n_customers
        cost_demands = [data.problem.costs[:,j] .* data.problem.demands[j] for j in 1:J]
        sorted_indices = [sortperm(cost_demands[j]) for j in 1:J]
        sorted_cost_demands = [cost_demands[j][sorted_indices[j]] for j in 1:J]
        
        obj_values = Vector{Float64}(undef, J)

        new(oracle_param, sorted_cost_demands, sorted_indices, J, obj_values)
    end

    UFLKnapsackOracle() = new()
end

function generate_cuts(oracle::UFLKnapsackOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol = 1e-9, time_limit = 3600.0, no_good_cut = false)
    tic = time()

    if no_good_cut && !all(x -> isapprox(0.5, abs(x-0.5), atol = tol), x_value)
        throw(AlgorithmException("generate_cuts: attempted to generate a no-good cut for a fractional solution"))
    end

    critical_facility = Vector{Int}(undef, oracle.J)
    for j in 1:oracle.J
        sorted_indices = oracle.sorted_indices[j]
        c_sorted = oracle.sorted_cost_demands[j]
        x_sorted = x_value[sorted_indices]

        # Find critical item and calculate contribution
        k = find_critical_item(c_sorted, x_sorted)

        # Calculate objective value contribution
        oracle.obj_values[j] = c_sorted[k] - (k > 1 ? sum((c_sorted[k] - c_sorted[i]) * x_sorted[i] for i in 1:k-1) : 0)

        if oracle.obj_values[j] >= t_value[j] * (1 + tol)
            critical_facility[j] = k
        else
            critical_facility[j] = oracle.oracle_param.add_only_violated_cuts ? -1 : k
        end
    end
    
    if get_sec_remaining(tic, time_limit) <= 0.0
        throw(TimeLimitException("Time limit reached during cut generation"))
    end

    customers = findall(x -> x != -1, critical_facility)

    # is_in_L should be determined by the sum of t's, must not individually
    is_in_L = sum(oracle.obj_values) >= sum(t_value) * (1 + tol) ? false : true

    hyperplanes = Vector{Hyperplane}()
    if !no_good_cut
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
    else
        h = Hyperplane(length(x_value), oracle.J)
        h.a_t = -ones(oracle.J)
        h.a_0 = sum(oracle.obj_values)
        for i=1:length(x_value)
            if x_value[i] > 0.5
                h.a_x[i] = -sum(oracle.obj_values) 
            end
        end
        push!(hyperplanes, h)
    end

    if is_in_L
        return !(oracle.oracle_param.slim) ? (true, hyperplanes, deepcopy(t_value)) : (true, [aggregate(hyperplanes)], deepcopy(t_value))
    end
    return !(oracle.oracle_param.slim) ? (false, hyperplanes, deepcopy(oracle.obj_values)) : (false, [aggregate(hyperplanes)], deepcopy(oracle.obj_values))
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