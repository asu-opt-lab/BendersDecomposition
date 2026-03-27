using Base.Threads: @threads
using OSQP, CPLEX, Gurobi

struct FacilityKnapsackInfo
    costs::Matrix{Float64}
    demands::Vector{Float64}
    capacity::Vector{Float64}
end

"""
    CFLKnapsackOracle

Oracle for generating Benders cuts for the capacitated facility location subproblem.

`CFLKnapsackOracle` solves the subproblem associated with a fixed first-stage
solution `x` and returns either an optimality cut or a feasibility cut. The
oracle supports generalized bound constraints, optional manual warm starts, and
a two-stage solve strategy in which the subproblem is first solved with a fast
configuration and then re-solved with a more reliable method when the objective
value is close to the current threshold.

When manual warm starts are enabled, the oracle uses a direct Gurobi model and
stores the previous primal and dual solution information so that it can be reused
across consecutive subproblem solves.
"""

mutable struct CFLKnapsackOracleParam <: AbstractOracleParam
    atol::Float64
    rtol::Float64

    function CFLKnapsackOracleParam(;atol = 1e-9, rtol = 1e-9)
        new(atol, rtol)
    end
end

"""
    CFLKnapsackOracle(data;
                      scen_idx=-1,
                      solver_param=Dict(...),
                      oracle_param=CFLKnapsackOracleParam(),
                      gbc=true,
                      manual_warm=false,
                      prev_y=zeros(...),
                      prev_dual=zeros(...))

Construct a cut-generation oracle for the capacitated facility location subproblem.

# Keyword Arguments
- `solver_param`: Solver attributes used when the model is created through the
  generic JuMP interface.
- `oracle_param`: Numerical tolerances used by the oracle when deciding whether
  a cut is violated.
- `gbc`: If `true`, generalized bound constraints are enforced by tightening the
  upper bounds of the assignment variables using the current `x` value.
- `manual_warm`: If `true`, the oracle uses a direct Gurobi model and manually
  loads primal and dual warm-start information.
- `prev_y`: Primal warm-start values for the assignment variables.
- `prev_dual`: Dual warm-start values for the stored constraints.

# Notes
- A direct Gurobi model is required when manual warm starts are enabled.
"""

mutable struct CFLKnapsackOracle <: AbstractTypicalOracle
    oracle_param::CFLKnapsackOracleParam
    solver_param::Dict{String,Any}

    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    facility_knapsack_info::FacilityKnapsackInfo
    gbc::Bool
    manual_warm::Bool
    prev_y::Vector{Float64}
    prev_dual::Vector{Float64}

    function CFLKnapsackOracle(data::Data; 
                               scen_idx=-1, 
                               solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9),
                               oracle_param::CFLKnapsackOracleParam = CFLKnapsackOracleParam(),
                               gbc::Bool = true,
                               manual_warm::Bool = false,
                               prev_y = zeros(data.problem.n_facilities*data.problem.n_customers),
                               prev_dual = zeros(data.problem.n_customers + data.problem.n_facilities*2))
        @debug "Building classical oracle"

        # A direct model is required when primal/dual warm starts are set manually.
        if manual_warm
            optimizer = optimizer_with_attributes(Gurobi.Optimizer, "Method" => 6, "PDHGGPU" => 1, "Crossover" => 1, "Threads" => 7, "LPWarmStart" => 0)
            model = direct_model(optimizer)
        else
            model = Model()
            assign_attributes!(model, solver_param)
        end

        # Define coupling variables and constraints
        @variable(model, 0 <= x[1:data.dim_x] <= 1)
        @constraint(model, fix_x, x .== 0)

        facility_knapsack_info = scen_idx == -1 ? FacilityKnapsackInfo(data.problem.costs, data.problem.demands, data.problem.capacities) : FacilityKnapsackInfo(data.problem.costs, data.problem.demands[scen_idx], data.problem.capacities)
        
        new(oracle_param, solver_param, model, fix_x, facility_knapsack_info, gbc, manual_warm, prev_y, prev_dual)
    end
    
    CFLKnapsackOracle() = new()
end

function generate_cuts(oracle::CFLKnapsackOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    # The subproblem is first solved with the default inexact configuration.
    sub_exact_solved = false
    s_time = time()
    set_time_limit_sec(oracle.model, time_limit)
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)

    # Apply generalized bound constraints by tightening y upper bounds to x.
    if oracle.gbc
        I, J = size(oracle.facility_knapsack_info.costs)
        for i in 1:I, j in 1:J
            set_upper_bound(oracle.model[:y][i, j], x_value[i])
        end
    end

    # Manually load primal and dual warm-start information.
    if oracle.manual_warm
        constrs = vcat(
            JuMP.all_constraints(oracle.model, JuMP.AffExpr, MOI.EqualTo{Float64}),
            JuMP.all_constraints(oracle.model, JuMP.AffExpr, MOI.LessThan{Float64}),
            JuMP.all_constraints(oracle.model, JuMP.AffExpr, MOI.GreaterThan{Float64}),
        )

        grb = JuMP.unsafe_backend(oracle.model)
        Gurobi.GRBupdatemodel(grb) 

        vars = all_variables(oracle.model)
        for (v, p_val) in zip(vars, vcat(x_value, oracle.prev_y))
            MOI.set(oracle.model, Gurobi.VariableAttribute("PStart"), v, p_val)
        end

        for (c, d_val) in zip(constrs, oracle.prev_dual)
            MOI.set(oracle.model, Gurobi.ConstraintAttribute("DStart"), c, d_val)
        end
    end

    optimize!(oracle.model)

    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end

    # If the inexact solve is close to the current threshold, re-solve with
    # Gurobi's automatic method selection to obtain a more reliable value.
    sub_obj_val = objective_value(oracle.model)
    if sub_obj_val < t_value[1] * (1 + oracle.oracle_param.rtol / tol_normalize)
        set_optimizer_attribute(oracle.model, "Method", -1)
        optimize!(oracle.model)
        sub_exact_solved = true
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

        # Store primal/dual information for the next warm start.
        oracle.manual_warm && (oracle.prev_y = vec(value.(oracle.model[:y])); oracle.prev_dual = vcat(dual.(oracle.model[:fix_x]), μ, dual.(oracle.model[:capacity])))

        # Restore the original PDHG setting after the exact re-solve.
        sub_exact_solved && set_optimizer_attribute(oracle.model, "Method", 6)
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