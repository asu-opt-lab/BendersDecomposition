using Base.Threads: @threads
using OSQP, CPLEX

struct FacilityKnapsackInfo
    costs::Matrix{Float64}
    demands::Vector{Float64}
    capacity::Vector{Float64}
end

mutable struct CFLKnapsackOracleParam <: AbstractOracleParam
    atol::Float64
    rtol::Float64

    function CFLKnapsackOracleParam(;atol = 1e-9, rtol = 1e-9)
        new(atol, rtol)
    end
end

mutable struct CFLKnapsackOracle <: AbstractTypicalOracle
    oracle_param::CFLKnapsackOracleParam
    solver_param::Dict{String,Any}

    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    facility_knapsack_info::FacilityKnapsackInfo

    function CFLKnapsackOracle(data::Data; 
                               scen_idx=-1, 
                               solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9),
                               oracle_param::CFLKnapsackOracleParam = CFLKnapsackOracleParam())
        @debug "Building classical oracle"
        model = Model()

        # Define coupling variables and constraints
        @variable(model, 0 <= x[1:data.dim_x] <= 1)
        @constraint(model, fix_x, x .== 0)

        facility_knapsack_info = scen_idx == -1 ? FacilityKnapsackInfo(data.problem.costs, data.problem.demands, data.problem.capacities) : FacilityKnapsackInfo(data.problem.costs, data.problem.demands[scen_idx], data.problem.capacities)

        assign_attributes!(model, solver_param)
        
        new(oracle_param, solver_param, model, fix_x, facility_knapsack_info)
    end
    
    CFLKnapsackOracle() = new()
end

function generate_cuts(oracle::CFLKnapsackOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    s_time = time()
    set_time_limit_sec(oracle.model, time_limit)
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    optimize!(oracle.model)

    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end

    status = dual_status(oracle.model)

    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(oracle.model)
        μ = dual.(oracle.model[:demand])
        a_t = [-1.0] 
        
        # Get facility knapsack info
        costs = oracle.facility_knapsack_info.costs
        demands = oracle.facility_knapsack_info.demands
        capacity = oracle.facility_knapsack_info.capacity

        # Calculate KP values for each facility
        KP_values = Vector{Float64}(undef, length(capacity))
        for i in 1:length(capacity)
            KP_values[i] = calculate_KP_value(costs[i,:], demands, capacity[i], μ)
        end

        a_x = KP_values # Vector{Float64}
        a_0 = sum(μ) 
        println("Exact time: $(time() - s_time)")
        if sub_obj_val >= t_value[1] * (1 + oracle.oracle_param.rtol / tol_normalize)
            return false, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        else
            return true, [Hyperplane(a_x, a_t, a_0)], deepcopy(t_value)
        end
        
    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(oracle.model)
            dual_sub_obj_val = dual_objective_value(oracle.model)
            @info "dual_sub_obj_val = $dual_sub_obj_val"
            
            a_x = dual.(oracle.fixed_x_constraints)
            a_t = [0.0]
            a_0 = dual_sub_obj_val - a_x' * x_value 
            if dual_sub_obj_val >= oracle.oracle_param.atol / tol_normalize
                return false, [Hyperplane(a_x, a_t, a_0)], [Inf]
            else
                return true, [Hyperplane(a_x, a_t, a_0)], [Inf]
            end
        end
        
    else
        println(status)
        throw(UnexpectedModelStatusException("ClassicalOracle: $(status)"))
    end
end

function calculate_KP_value(costs::Vector{Float64}, demands::Vector{Float64}, capacity::Float64, μ::Vector{Float64})
    n = length(demands)
    
    # ratios = Vector{Tuple{Int,Float64}}(undef, n)
    ratios = [(i, (costs[i] * demands[i] - μ[i]) / demands[i]) for i in 1:n if (costs[i] * demands[i] - μ[i]) < 0]
    
    sort!(ratios, by=x->x[2])
    
    kp_value = 0.0
    remaining_capacity = capacity
    z = zeros(n)

    for (i, _) in ratios
        if remaining_capacity >= demands[i]
            kp_value += costs[i] * demands[i] - μ[i]
            remaining_capacity -= demands[i]
            z[i] = 1.0
        else
            fraction = remaining_capacity / demands[i]
            kp_value += (costs[i]*demands[i] - μ[i]) * fraction
            z[i] = fraction
            break
        end
    end

    return kp_value
end

# no multi pass improvment
function UB_approximation(data, x; rtol =1e-6, tol=1e-6)
    costs, demands, capacities = data.problem.costs, data.problem.demands, data.problem.capacities

    m, n = size(costs)
    total_demands, usuable_capa = sum(demands), capacities .* x

    # Uniform distribution
    s = map(i -> min(x[i], usuable_capa[i] / total_demands), 1:m); S = sum(s)
    @assert S ≥ 1 - tol "Uniform distribution cannot be applied (S<1)"
    w = s ./ S # scaling to satisfy Σ_j y_ij = 1 ∀ j
    y = repeat(w, 1, n) # distriubte w_i to each row

    # Improve obj
    residual = usuable_capa .- total_demands .* w
    x_ascen_cost_idxes = [sortperm(costs[:, j]) for j in 1:n]
    for j in 1:n
        d = demands[j]; (d ≤ tol && continue)

        idxes = x_ascen_cost_idxes[j] # ascending idxes of x w.r.t. costs
        receiver = 1; donor = m # chepeast, most expensive idx

        while receiver < donor
            i_receive = idxes[receiver]; i_donate = idxes[donor] # corresponding facilities 

            c_receive = costs[i_receive, j]; c_donate = costs[i_donate, j]
            c_donate > c_receive || break # if c_receive is large than no improvment at current jth col

            # Compute cap of receiver
            arc_slack_receive = max(0.0, x[i_receive] - y[i_receive, j]) # arc cap
            cap_slack_receive = max(0.0, residual[i_receive] / d) # capa cap
            receive = min(arc_slack_receive, cap_slack_receive)
            receive ≤ tol && (receiver += 1; continue) # if no remaining choose other receiver

            # Compute cap of donor
            donate = y[i_donate, j]
            donate ≤ tol && (donor -= 1; continue) # if cannot donate, choose other donor

            # if receive = donate = 0
            Δ = min(receive, donate)
            if Δ ≤ tol
                receive ≤ tol && (receiver += 1)
                donate ≤ tol && (donor -= 1)
                continue
            end

            # Update y, residual
            y[i_donate, j] -= Δ; y[i_receive, j] += Δ
            residual[i_receive] -= d * Δ; residual[i_donate] += d * Δ

            # Check y = 0
            y[i_donate, j] ≤ tol && (donor -= 1)

            # Update cap of donor
            arc_slack_receive = max(0.0, x[i_receive] - y[i_receive, j])
            cap_slack_receive = max(0.0, residual[i_receive] / d)
            receive = min(arc_slack_receive, cap_slack_receive)
            receive ≤ tol && (receiver += 1)
        end
    end

    @assert all(abs.(sum(y, dims=1)[:] .- 1.0) .≤ tol)
    @assert all(y .≥ -tol .&& y .≤ x .+ tol)
    @assert all(i -> sum(demands .* y[i, :]) ≤ capacities[i]*x[i] + tol + rtol*capacities[i]*x[i], 1:length(x))

    # return sum((costs.* demands') .* y) * 2
    return sum((costs.* demands') .* y)
end