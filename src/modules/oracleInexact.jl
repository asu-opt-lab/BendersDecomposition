export InexactOracle

"""
    InexactOracle

Oracle for generating cuts by solving subproblems either approximately or exactly.

`InexactOracle` is designed for use with solvers such as CuOpt or OSQP. It first
solves the subproblem with an inexact solver and, when necessary, switches to an
exact solver to validate the result.

Because CuOpt does not support `JuMP.dual`, the oracle relies on an equivalent
dual formulation when dual information is required.
"""

struct FacilityKnapsackInfo
    costs::Matrix{Float64}
    demands::Vector{Float64}
    capacity::Vector{Float64}
end

"""
    InexactOracle(data, exact_solver_param, inexact_solver_param;
                  scen_idx=-1,
                  oracle_param=InexactOracleParam(),
                  exact_solve=false)

Construct an oracle that solves subproblems inexactly by default and optionally
switches to an exact solver when higher accuracy is needed.

# Arguments
- `data`: Problem data.
- `exact_solver_param`: Solver attributes for the exact solver.
- `inexact_solver_param`: Solver attributes for the inexact solver.

# Keyword Arguments
- `scen_idx`: Scenario index. If `-1`, the full demand vector is used.
- `oracle_param`: Tolerance parameters for the oracle.
- `exact_solve`: Initial flag indicating whether the exact solver is active.

# Notes
- If the inexact solver is not CuOpt, fixing constraints for `x` are explicitly
  created and stored in `fixed_x_constraints`.
- If the inexact solver is CuOpt, those fixing constraints are omitted and the
  variable fixing is handled through the model objective structure instead.
"""

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

        # If the solver is not CuOpt, explicitly maintain fixing constraints for x.
        # CuOpt uses a different mechanism, so these constraints are omitted.
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

"""
    generate_cuts(oracle::InexactOracle, x_value, t_value;
                  tol_normalize=1.0,
                  time_limit=3600)

Generate cuts for the current master solution by solving the oracle subproblem.

The subproblem is first solved with the inexact solver stored in `oracle`. If the
result is close enough to the current `t_value`, the model is re-solved with the
exact solver to avoid accepting an inaccurate cut.

# Notes
- For CuOpt, dual quantities are recovered from the dual reformulation because
  `JuMP.dual` is not supported.
- For non-CuOpt solvers, fixed-`x` constraints are updated directly through their
  right-hand sides.
- If the subproblem reaches the time limit, a `TimeLimitException` is thrown.
- If the solver returns an unexpected status, an `UnexpectedModelStatusException`
  is thrown.
"""

function generate_cuts(oracle::InexactOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    model = oracle.model; oracle_param = oracle.oracle_param
    solver = MOI.get(model, MOI.SolverName())

    # The way x is injected into the subproblem depends on whether the solver is CuOpt.
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

    # If the inexact objective is sufficiently close to improving the current t-value,
    # re-solve the subproblem with the exact solver.
    if sub_obj_val < t_value[1] * (1 + oracle_param.rtol / tol_normalize)
        assign_attributes!(model, oracle.exact_solver_param)
        optimize!(model)
        sub_obj_val = objective_value(model)
        oracle.exact_solve = true
    end

    # primal status should be checked for cuOpt
    status = (solver == "cuOpt") ? primal_status(model) : dual_status(model)

    if status == FEASIBLE_POINT
        # μ is equivalent to p of dual model (see model.jl of cflp)
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
        
        # Reassign inexact solver
        oracle.exact_solve && assign_attributes!(model, oracle.inexact_solver_param)
        if sub_obj_val >= t_value[1] * (1 + oracle_param.rtol / tol_normalize)
            # Obj of subproblem can be returned if subproblem is solved exactly
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
            
            a_x = dual.(oracle.fixed_x_constraints)
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