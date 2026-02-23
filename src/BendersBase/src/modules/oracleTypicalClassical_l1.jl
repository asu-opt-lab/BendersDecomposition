export ClassicalOracle

using JuMP
using MathOptInterface
const MOI = MathOptInterface

# ============================================================================
# L1 feasibility normalization state (no GBC, no master-copy machinery)
# ============================================================================

mutable struct L1State
    # L1 normalization variables (persistent)
    z0::VariableRef
    mu_plus::Vector{VariableRef}
    mu_minus::Vector{VariableRef}

    # Constraints involving x that will be relaxed in L1 mode
    l1_geq_constraints::Vector{ConstraintRef}   # >= : add +z0
    l1_leq_constraints::Vector{ConstraintRef}   # <= : add -z0
    l1_eq_constraints::Vector{ConstraintRef}    # == : add +mu_plus - mu_minus

    # Restore objective after L1 mode
    original_objective::Any
    original_objective_sense::MOI.OptimizationSense

    # z0 = sum(mu_plus) + sum(mu_minus) (optional)
    z0_sum_constraint::Union{ConstraintRef, Nothing}
end

mutable struct ClassicalOracle <: AbstractTypicalOracle
    oracle_param::BasicOracleParam

    model::Model
    fixed_x_constraints::Vector{ConstraintRef}

    # L1 feasibility state
    l1_state::L1State

    # ------------------------------------------------------------------------
    # Auto-decompose constructor (model already built elsewhere)
    # ------------------------------------------------------------------------
    function ClassicalOracle(oracle_param::BasicOracleParam,
                             model::Model,
                             fixed_x_constraints::Vector{ConstraintRef})

        # Recover x variables from fixing constraints (expects each is x_i == rhs)
        x_vars = VariableRef[]
        for con in fixed_x_constraints
            func = constraint_object(con).func

            if func isa VariableRef
                push!(x_vars, func)

            elseif func isa AffExpr
                # Accept the single-variable form: (1.0)*x + 0 == rhs
                if length(func.terms) == 1 && iszero(func.constant)
                    (v, coef) = first(func.terms)
                    if coef == 1.0
                        push!(x_vars, v)
                    else
                        throw(ArgumentError("Fixing constraint has coefficient $coef; expected 1.0"))
                    end
                else
                    throw(ArgumentError("Fixing constraint is not a single-variable equality (AffExpr has multiple terms or nonzero constant)."))
                end
            else
                throw(ArgumentError("Unsupported fixing constraint function type: $(typeof(func))."))
            end
        end

        l1_state = _setup_l1_state(model, x_vars, fixed_x_constraints)
        new(oracle_param, model, fixed_x_constraints, l1_state)
    end

    # ------------------------------------------------------------------------
    # First-code-style constructor (build model with x and fix_x, user customizes)
    # ------------------------------------------------------------------------
    function ClassicalOracle(data::Data;
                             scen_idx::Int = -1,
                             solver_param::Dict{String,Any} = Dict(
                                 "solver" => "CPLEX",
                                 "CPX_PARAM_EPRHS" => 1e-9,
                                 "CPX_PARAM_NUMERICALEMPHASIS" => 1,
                                 "CPX_PARAM_EPOPT" => 1e-9
                             ),
                             oracle_param::BasicOracleParam = BasicOracleParam(),
                             build_submodel!::Function = (model, data, scen_idx, x) -> nothing)

        @debug "Building classical oracle (simple + L1 feasibility)"
        model = Model()

        # Coupling variables and fixing constraints (same style as your original)
        @variable(model, x[1:data.dim_x])
        @constraint(model, fix_x, x .== 0)

        assign_attributes!(model, solver_param)

        # User hook: add subproblem objective/constraints using x
        build_submodel!(model, data, scen_idx, x)

        l1_state = _setup_l1_state(model, x, fix_x)
        new(oracle_param, model, fix_x, l1_state)
    end

    ClassicalOracle() = new()
end

# ============================================================================
# L1 setup helpers
# ============================================================================

"""
    _constraint_involves_x(con, x_indices) -> Bool

Return true if constraint involves any x variable.
"""
function _constraint_involves_x(con::ConstraintRef, x_indices::Set)
    func = constraint_object(con).func

    if func isa VariableRef
        return func.index in x_indices
    elseif func isa AffExpr
        for (var, _) in func.terms
            if var.index in x_indices
                return true
            end
        end
    end
    return false
end

"""
    _setup_l1_state(model, x_vars, fix_x_constraints) -> L1State

Scan constraints involving x (excluding fix_x constraints) and prepare z0 and mu variables.
"""
function _setup_l1_state(model::Model,
                         x_vars::Vector{VariableRef},
                         fix_x_constraints::Vector{ConstraintRef})

    x_indices = Set(v.index for v in x_vars)
    fix_x_indices = Set(c.index for c in fix_x_constraints)

    # Save objective for restoration
    original_objective = objective_function(model)
    original_objective_sense = objective_sense(model)

    l1_geq_constraints = ConstraintRef[]
    l1_leq_constraints = ConstraintRef[]
    l1_eq_constraints  = ConstraintRef[]

    # Identify constraints that involve x (skip the fixing constraints)
    for con in all_constraints(model, include_variable_in_set_constraints=false)
        if con.index in fix_x_indices
            continue
        end

        if !_constraint_involves_x(con, x_indices)
            continue
        end

        set = constraint_object(con).set
        if set isa MOI.GreaterThan
            push!(l1_geq_constraints, con)
        elseif set isa MOI.LessThan
            push!(l1_leq_constraints, con)
        elseif set isa MOI.EqualTo
            push!(l1_eq_constraints, con)
        end
    end

    # z0 >= 0
    z0 = @variable(model, lower_bound = 0)
    set_name(z0, "z0")

    # mu+/mu- for each equality constraint
    mu_plus = VariableRef[]
    mu_minus = VariableRef[]
    n_eq = length(l1_eq_constraints)

    z0_sum_constraint = nothing
    if n_eq > 0
        for i in 1:n_eq
            mu_p = @variable(model, lower_bound = 0)
            mu_m = @variable(model, lower_bound = 0)
            set_name(mu_p, "mu_plus[$i]")
            set_name(mu_m, "mu_minus[$i]")
            push!(mu_plus, mu_p)
            push!(mu_minus, mu_m)
        end

        z0_sum_constraint = @constraint(model, z0 == sum(mu_plus) + sum(mu_minus))
        set_name(z0_sum_constraint, "z0_sum")
    end

    @debug "L1 setup complete" n_geq=length(l1_geq_constraints) n_leq=length(l1_leq_constraints) n_eq=n_eq

    return L1State(
        z0,
        mu_plus,
        mu_minus,
        l1_geq_constraints,
        l1_leq_constraints,
        l1_eq_constraints,
        original_objective,
        original_objective_sense,
        z0_sum_constraint
    )
end

# ============================================================================
# L1 mode toggling
# ============================================================================

"""
    _enter_l1_mode!(oracle)

Activate z0 / mu coefficients and set objective to Min z0.
Requires `set_normalized_coefficient(con, var, val)` to exist in your codebase.
"""
function _enter_l1_mode!(oracle::ClassicalOracle)
    s = oracle.l1_state

    # >= constraints: +z0
    for con in s.l1_geq_constraints
        set_normalized_coefficient(con, s.z0, 1.0)
    end

    # <= constraints: -z0
    for con in s.l1_leq_constraints
        set_normalized_coefficient(con, s.z0, -1.0)
    end

    # == constraints: +mu_plus - mu_minus
    for (i, con) in enumerate(s.l1_eq_constraints)
        set_normalized_coefficient(con, s.mu_plus[i], 1.0)
        set_normalized_coefficient(con, s.mu_minus[i], -1.0)
    end

    @objective(oracle.model, Min, s.z0)
end

"""
    _exit_l1_mode!(oracle)

Deactivate z0 / mu coefficients and restore original objective.
"""
function _exit_l1_mode!(oracle::ClassicalOracle)
    s = oracle.l1_state

    for con in s.l1_geq_constraints
        set_normalized_coefficient(con, s.z0, 0.0)
    end
    for con in s.l1_leq_constraints
        set_normalized_coefficient(con, s.z0, 0.0)
    end
    for (i, con) in enumerate(s.l1_eq_constraints)
        set_normalized_coefficient(con, s.mu_plus[i], 0.0)
        set_normalized_coefficient(con, s.mu_minus[i], 0.0)
    end

    set_objective_function(oracle.model, s.original_objective)
    set_objective_sense(oracle.model, s.original_objective_sense)
end

# ============================================================================
# L1 feasibility cut generation
# ============================================================================

function _generate_l1_feasibility_cut(oracle::ClassicalOracle,
                                      x_value::Vector{Float64},
                                      tol_normalize::Float64,
                                      time_limit::Float64)

    _enter_l1_mode!(oracle)

    set_time_limit_sec(oracle.model, time_limit)

    # Keep linking constraints x == x̂
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)

    optimize!(oracle.model)

    if termination_status(oracle.model) == TIME_LIMIT
        _exit_l1_mode!(oracle)
        throw(TimeLimitException("Time limit reached during L1 feasibility cut generation"))
    end

    if termination_status(oracle.model) != OPTIMAL
        _exit_l1_mode!(oracle)
        throw(UnexpectedModelStatusException(
            "ClassicalOracle L1 mode returned $(termination_status(oracle.model)); expected OPTIMAL."
        ))
    end

    z0_val = objective_value(oracle.model)

    a_x = dual.(oracle.fixed_x_constraints)
    a_t = [0.0]
    a_0 = z0_val - a_x' * x_value

    _exit_l1_mode!(oracle)

    if z0_val >= oracle.oracle_param.zero_tol / tol_normalize
        return false, [Hyperplane(a_x, a_t, a_0)], [Inf]
    else
        return true,  [Hyperplane(a_x, a_t, a_0)], [Inf]
    end
end

# ============================================================================
# Main cut generation (optimality if feasible; otherwise L1 feasibility)
# ============================================================================

function generate_cuts(oracle::ClassicalOracle,
                       x_value::Vector{Float64},
                       t_value::Vector{Float64};
                       tol_normalize = 1.0,
                       time_limit = 3600)

    set_time_limit_sec(oracle.model, time_limit)
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)

    optimize!(oracle.model)

    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end

    status = dual_status(oracle.model)

    # --------------------------
    # Optimality cut
    # --------------------------
    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(oracle.model)

        a_x = dual.(oracle.fixed_x_constraints)
        a_t = [-1.0]
        a_0 = sub_obj_val - a_x' * x_value

        if sub_obj_val >= t_value[1] * (1 + oracle.oracle_param.rtol) + oracle.oracle_param.atol / tol_normalize
            return false, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        else
            return true,  [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        end
    elseif status == INFEASIBILITY_CERTIFICATE || termination_status(oracle.model) == INFEASIBLE
        # Use L1 normalization for feasibility cuts
        return _generate_l1_feasibility_cut(oracle, x_value, tol_normalize, time_limit)
    else
        throw(UnexpectedModelStatusException("ClassicalOracle: $(status). This is likely a numerical issue. Please try using other oracles, such as unified oracle or pareto oracle."))
    end
end