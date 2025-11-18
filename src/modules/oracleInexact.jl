export InexactOracle

struct FacilityKnapsackInfo
    costs::Matrix{Float64}
    demands::Vector{Float64}
    capacity::Vector{Float64}
end

mutable struct InexactOracleParam <: AbstractOracleParam
    atol::Float64
    rtol::Float64

    function InexactOracleParam(;atol = 1e-9, rtol = 1e-9)
        new(atol, rtol)
    end
end

mutable struct InexactOracle <: AbstractTypicalOracle
    data::Data
    oracle_param::InexactOracleParam
    exact_solver_param::Dict{String,Any}
    inexact_solver_param::Dict{String,Any}

    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    facility_knapsack_info::FacilityKnapsackInfo

    exact_solve::Bool

    function InexactOracle(data, exact_solver_param, inexact_solver_param; scen_idx=-1, oracle_param::InexactOracleParam = InexactOracleParam(), exact_solve = false)
        model = Model()

        if inexact_solver_param["solver"] != "cuOpt"
            @variable(model, 0 <= x[1:data.dim_x] <= 1)
            @constraint(model, fix_x, x .== 0)
        else
            fix_x = JuMP.ConstraintRef[]
        end

        facility_knapsack_info = scen_idx == -1 ? FacilityKnapsackInfo(data.problem.costs, data.problem.demands, data.problem.capacities) : FacilityKnapsackInfo(data.problem.costs, data.problem.demands[scen_idx], data.problem.capacities)

        assign_attributes!(model, inexact_solver_param)
        
        new(data, oracle_param, exact_solver_param, inexact_solver_param, model, fix_x, facility_knapsack_info, exact_solve)
    end

    InexactOracle() = new()
end

function generate_cuts(oracle::InexactOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    model = oracle.model; oracle_param = oracle.oracle_param
    solver = MOI.get(model, MOI.SolverName())

    if solver == "cuOpt"
        set_objective_coefficient.(model, model[:s], x_value)
    else
        set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    end

    set_time_limit_sec(model, time_limit)
    optimize!(model)

    if termination_status(model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end
    sub_obj_val = objective_value(model)

    # dual feasibility
    data = oracle.data
    I, J = data.problem.n_facilities, data.problem.n_customers
    cost_demands = data.problem.costs .* data.problem.demands'
    p = value.(model[:p]); l = value.(model[:l]); d = value.(model[:d]); s = value.(model[:s])
    cons1 = [isapprox(sum(d[i,:]) + l[i]*data.problem.capacities[i] + s[i], 0.0; atol=1e-9) for i in 1:I]
    cons2 = [(p[j] - d[i,j] - l[i] * data.problem.demands[j] <= cost_demands[i,j] + 1e-9)  for i in 1:I, j in 1:J]
    println("# of cons1 violation: count(cons1), # of cons2 violation: count(cons2)")

    if sub_obj_val < t_value[1] * (1 + oracle_param.rtol / tol_normalize)
        println("obj of solver: $(objective_value(model)), t_value: $(t_value[1])")
        assign_attributes!(model, oracle.exact_solver_param)
        optimize!(model)
        sub_obj_val = objective_value(model)
        oracle.exact_solve = true
    end

    status = (solver == "cuOpt") ? primal_status(model) : dual_status(model)

    if status == FEASIBLE_POINT
        μ = (solver == "cuOpt") ? value.(model[:p]) : dual.(model[:demand])
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
        
        oracle.exact_solve && assign_attributes!(model, oracle.inexact_solver_param)
        if sub_obj_val >= t_value[1] * (1 + oracle_param.rtol / tol_normalize)
            oracle.exact_solve && (oracle.exact_solve = false; return false, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val])
            return false, [Hyperplane(a_x, a_t, a_0)], [NaN]
        else
            return true, [Hyperplane(a_x, a_t, a_0)], deepcopy(t_value)
        end
        
    elseif status == INFEASIBILITY_CERTIFICATE
        exist_dual = (solver == "cuOpt") ? has_values(model) : has_duals(model)
        if exist_dual
            dual_sub_obj_val = (solver == "cuOpt") ? sub_obj_val : dual_objective_value(model)
            @info "dual_sub_obj_val = $dual_sub_obj_val"
            
            a_x = dual.(typical_oracle.fixed_x_constraints)
            a_t = [0.0]
            a_0 = dual_sub_obj_val - a_x' * x_value 
            if dual_sub_obj_val >= oracle_param.atol / tol_normalize
                return false, [Hyperplane(a_x, a_t, a_0)], [Inf]
            else
                return true, [Hyperplane(a_x, a_t, a_0)], [Inf]
            end
        end
    else
        @info("termination_status: $(termination_status(model)), primal: $(primal_status(model)), dual: $status")
        throw(UnexpectedModelStatusException("InexactOracle: $(status)"))
    end
end