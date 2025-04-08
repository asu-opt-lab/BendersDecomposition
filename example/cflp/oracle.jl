struct FacilityKnapsackInfo
    costs::Matrix{Float64}
    demands::Vector{Float64}
    capacity::Vector{Float64}
end

mutable struct CFLKnapsackOracle <: AbstractTypicalOracle
    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    other_constraints::Vector{ConstraintRef}
    facility_knapsack_info::FacilityKnapsackInfo

    function CFLKnapsackOracle(data::Data; scen_idx=-1)
        @debug "Building classical oracle"
        model = Model()

        # Define coupling variables and constraints
        @variable(model, x[1:data.dim_x])
        @constraint(model, fix_x, x .== 0)

        other_constr = Vector{ConstraintRef}()

        facility_knapsack_info = scen_idx == -1 ? FacilityKnapsackInfo(data.problem.costs, data.problem.demands, data.problem.capacities) : FacilityKnapsackInfo(data.problem.costs, data.problem.demands[scen_idx], data.problem.capacities)

        new(model, fix_x, other_constr, facility_knapsack_info)
    end
    
    CFLKnapsackOracle() = new()
end


function generate_cuts(oracle::CFLKnapsackOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol=1e-6)
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    optimize!(oracle.model)
    status = dual_status(oracle.model)

    h = Hyperplane(length(x_value),length(t_value))

    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(oracle.model)

        if sub_obj_val >= t_value[1] + tol
            μ = dual.(oracle.model[:demand])
            h.a_t = [-1.0] 
            
            # Get facility knapsack info
            costs = oracle.facility_knapsack_info.costs
            demands = oracle.facility_knapsack_info.demands
            capacity = oracle.facility_knapsack_info.capacity

            # Calculate KP values for each facility
            KP_values = Vector{Float64}(undef, length(capacity))
            for i in 1:length(capacity)
                KP_values[i] = calculate_KP_value(costs[i,:], demands, capacity[i], μ)
            end

            h.a_x = KP_values # Vector{Float64}
            h.a_0 = sum(μ) 
            return false, [h], [sub_obj_val]
        else
            return true, [h], t_value
        end
        
    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(oracle.model)
            h.a_x = dual.(oracle.fixed_x_constraints)
            h.a_t = [0.0]
            h.a_0 = dual.(oracle.other_constraints)'*normalized_rhs.(oracle.other_constraints)
            return false, [h], [Inf]
        else
            throw(ErrorException("CFLKnapsackOracle oracle: Infeasible subproblem has no dual solution"))
        end
        
    else
        throw(ErrorException("CFLKnapsackOracle oracle: Unexpected dual status $status"))
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
