export ClassicalOracle, ClassicalOracleParam

const ClassicalOracleParam = BasicOracleParam

mutable struct L1State
    # L1 normalization fields (persistent)
    z0::VariableRef                                    # L1 penalty variable
    mu_plus::Vector{VariableRef}                       # Positive deviation for == constraints
    mu_minus::Vector{VariableRef}                      # Negative deviation for == constraints

    # Constraint tracking for L1 mode
    l1_geq_constraints::Vector{ConstraintRef}          # >= constraints involving x (add +z0)
    l1_leq_constraints::Vector{ConstraintRef}          # <= constraints involving x (add -z0)
    l1_eq_constraints::Vector{ConstraintRef}           # == constraints involving x (use mu)

    # State for restoration
    original_objective::Union{AffExpr, VariableRef, Real}
    original_objective_sense::MOI.OptimizationSense

    # Auxiliary constraint for z0 = sum(mu+ + mu-)
    z0_sum_constraint::Union{ConstraintRef, Nothing}

    # GBC base bounds (for L1 feasibility mode)
    gbc_base_lower::Vector{Float64}
    gbc_base_upper::Vector{Float64}
end

mutable struct ClassicalOracle <: AbstractTypicalOracle
    
    param::ClassicalOracleParam

    model::Model
    fixed_x_constraints::Vector{ConstraintRef}

    gbc_lhs::Vector{VariableRef}
    gbc_rhs::Vector{Union{VariableRef, AffExpr}}
    gbc_sense::Vector{GBCBoundType}
    l1_state::L1State


    function ClassicalOracle(data::AbstractData, master::Master; 
                            customize = customize_sub_model!,
                            scen_idx::Int=0, 
                            param::ClassicalOracleParam = ClassicalOracleParam())
    
            @debug "Building classical oracle"
            model = Model()

            # Copy the master's coupling variables into the submodel (with identical axes and symbols)
            x_copy = copy_variables!(model, master.x_tuple)

            # Collect all copied master variables and add linking constraint
            x = var_from_tuple(x_copy)
            @constraint(model, fix_x, x .== 0)

            # Build the submodel using user-defined customization, passing the copied variables
            result = customize(model, data, scen_idx; x_copy...)
            
            # Validate that all constraints are supported types (LP)
            _validate_constraint_types(model)
            
            # Parse the result to extract GBC information
            gbc_lhs, gbc_rhs, gbc_sense = _parse_gbc_result(result, x)

            # Capture base bounds for GBC variables (before any GBC bound updates)
            gbc_base_lower = Float64[has_lower_bound(var) ? lower_bound(var) : -Inf for var in gbc_lhs]
            gbc_base_upper = Float64[has_upper_bound(var) ? upper_bound(var) : Inf for var in gbc_lhs]

            # Create oracle instance with basic fields
            oracle = new()
            oracle.param = param
            oracle.model = model
            oracle.fixed_x_constraints = fix_x
            oracle.gbc_lhs = gbc_lhs
            oracle.gbc_rhs = gbc_rhs
            oracle.gbc_sense = gbc_sense
            
            # Setup L1 normalization state
            oracle.l1_state = _setup_l1_state(model, x, fix_x, gbc_base_lower, gbc_base_upper)
            
            return oracle
    end

    ClassicalOracle() = new()
end

"""
    _constraint_involves_x(con::ConstraintRef, x_indices::Set) -> Bool

Check if a constraint involves any of the master (x) variables.
Used to identify which constraints need L1 normalization.
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
    _setup_l1_state(model::Model, x_vars::Vector{VariableRef}, fix_x_constraints,
                    gbc_base_lower::Vector{Float64}, gbc_base_upper::Vector{Float64})

Setup persistent L1 normalization variables and identify constraints that involve x.
This is called once during construction.

Per L1.md formulation:
- For >= constraints: add z0 (coefficient +1.0 when active)
- For <= constraints: add z0 (coefficient -1.0 when active)
- For == constraints: use mu+/mu- pairs (B₀y + μ⁺ - μ⁻ = b₀ - A₀x̂)
- z0 = sum(mu+) + sum(mu-) (L1 norm relationship)
"""
function _setup_l1_state(model::Model, x_vars::Vector{VariableRef}, fix_x_constraints,
                         gbc_base_lower::Vector{Float64}, gbc_base_upper::Vector{Float64})
    # Build set of master variable indices and fixing constraint indices for fast lookup
    x_indices = Set(v.index for v in x_vars)
    fix_x_indices = Set(c.index for c in fix_x_constraints)

    # Store original objective (for restoration after L1 mode)
    original_objective = objective_function(model)
    original_objective_sense = objective_sense(model)

    # Initialize constraint tracking vectors
    l1_geq_constraints = ConstraintRef[]
    l1_leq_constraints = ConstraintRef[]
    l1_eq_constraints = ConstraintRef[]

    # Identify constraints involving x and categorize by type
    for con in all_constraints(model, include_variable_in_set_constraints=false)
        # Skip fixing constraints (they are for x = x*)
        if con.index in fix_x_indices
            continue
        end

        # Only process constraints that involve master variables
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

    # Create z0 penalty variable (z0 >= 0)
    z0 = @variable(model, lower_bound = 0)
    set_name(z0, "z0")

    # Create mu+/mu- pairs for each equality constraint
    mu_plus = VariableRef[]
    mu_minus = VariableRef[]
    n_eq = length(l1_eq_constraints)
    if n_eq > 0
        for i in 1:n_eq
            mu_p = @variable(model, lower_bound = 0)
            mu_m = @variable(model, lower_bound = 0)
            set_name(mu_p, "mu_plus[$i]")
            set_name(mu_m, "mu_minus[$i]")
            push!(mu_plus, mu_p)
            push!(mu_minus, mu_m)
        end

        # Create z0 sum constraint: z0 = sum(mu+) + sum(mu-)
        # This is always active, but only matters when in L1 mode
        z0_sum_constraint = @constraint(model, z0 == sum(mu_plus) + sum(mu_minus))
        set_name(z0_sum_constraint, "z0_sum")
    else
        z0_sum_constraint = nothing
    end

    @debug "L1 setup complete" n_geq=length(l1_geq_constraints) n_leq=length(l1_leq_constraints) n_eq=length(l1_eq_constraints)

    return L1State(
        z0,
        mu_plus,
        mu_minus,
        l1_geq_constraints,
        l1_leq_constraints,
        l1_eq_constraints,
        original_objective,
        original_objective_sense,
        z0_sum_constraint,
        gbc_base_lower,
        gbc_base_upper,
    )
end

"""
    _parse_gbc_result(result, x_vars) -> (gbc_lhs, gbc_rhs, gbc_sense)

Parse and validate the value returned by a user `customize` function and
convert it into the internal GBC representation.

# Returns
(gbc_lhs, gbc_rhs, gbc_sense)
- `gbc_lhs :: Vector{VariableRef}`  
  Left-hand-side variables for the generalized bound constraint. Every item
  must be a `VariableRef` that is **not** one of the copied master variables
  supplied in `x_vars`.
- `gbc_rhs :: Vector{Union{VariableRef, AffExpr}}`  
  Right-hand-side expressions. Each item must be either
  - a `VariableRef` drawn from `x_vars` (i.e., a *copied* master variable), or
  - an `AffExpr` expressed only in terms of variables from `x_vars` (and
    numeric constants).
- `gbc_sense :: Vector{GBCBoundType}`  
  A vector of bound senses for each constraint (for example `UpperBound`,
  `LowerBound`, or `Fixed`).

# Arguments
- `result` :: Any
  The raw return value from the user `customize` function. This function must
  parse and normalize `result` into the three-tuple described above.
- `x_vars` :: Vector{VariableRef}  
  A vector of submodel variables that are *copies* of the coupling master
  variables. Used only for validation.

# Validation
1. `gbc_lhs` must be a `Vector{VariableRef}`. None of the entries in `gbc_lhs`
   may be any of the `VariableRef`s listed in `x_vars`.
2. `gbc_rhs` must be a `Vector{Union{VariableRef, AffExpr}}`. Every `VariableRef`
   appearing in `gbc_rhs` must be one of the `VariableRef`s in `x_vars`.
3. `gbc_sense` must be a `Vector{GBCBoundType}` with the same length as `gbc_lhs`
   and `gbc_rhs`.

# Errors
If any validation rule is violated, raise a descriptive exception (for example
`ArgumentError`). Error messages should identify:
- which part (`gbc_lhs` / `gbc_rhs` / `gbc_sense`) failed validation, and
- the offending element(s) and why (e.g., "lhs contains copied master variable `x[3]`",
  or "rhs contains variable `y` not in `x_vars`").

# Examples
See the test suite for representative `result` shapes and the exact exception
types.
"""
function _parse_gbc_result(result, x_vars::Vector{VariableRef})
    # No GBC - must explicitly return nothing
    if result === nothing
        return VariableRef[], Union{VariableRef, AffExpr}[], GBCBoundType[]
    end
    
    # Format: (gbc_lhs, gbc_rhs, gbc_sense)
    if length(result) == 3
        gbc_lhs, gbc_rhs, gbc_sense = result
        
        # Validate lengths
        n = length(gbc_lhs)
        if length(gbc_rhs) != n || length(gbc_sense) != n
            throw(DimensionMismatch(
                "All GBC vectors must have the same length. " *
                "Got length(gbc_lhs) = $n, length(gbc_rhs) = $(length(gbc_rhs)), " *
                "length(gbc_sense) = $(length(gbc_sense))."
            ))
        end
        
        if isempty(gbc_lhs)
            @warn "GBC tuple returned but gbc_lhs is empty. No GBC constraints will be applied."
            return VariableRef[], Union{VariableRef, AffExpr}[], GBCBoundType[]
        end
        
        # Build set of valid master variable indices for fast lookup
        x_indices = Set(v.index for v in x_vars)
        
        # Validate gbc_lhs: must be a vector of VariableRef and must NOT contain copied master variables (should be subproblem variables)
        for (i, lhs_var) in enumerate(gbc_lhs)
            if !(lhs_var isa VariableRef)
                throw(ArgumentError("gbc_lhs[$i] is not a VariableRef (got $(typeof(lhs_var))). Each lhs entry must be a VariableRef."))
            end
            if lhs_var.index in x_indices
                throw(ArgumentError(
                    "gbc_lhs[$i] contains a copied master variable $lhs_var. " *
                    "Each item of gbc_lhs should only contain a single subproblem variable. " *
                    "The GBC relation should be: subproblem_var <=/>=/== affine_func(master_vars)."
                ))
            end
        end
        
        # Validate gbc_rhs: must ONLY contain copied master variables in affine form
        for (i, rhs) in enumerate(gbc_rhs)
            if rhs isa VariableRef
                if !(rhs.index in x_indices)
                    throw(ArgumentError(
                        "gbc_rhs[$i] contains a variable $(rhs) that is not a copied master variable. " *
                        "gbc_rhs should only reference copied master variables passed to customize(). "
                    ))
                end
            elseif rhs isa AffExpr  # AffExpr
                for (var, _) in rhs.terms
                    if !(var.index in x_indices)
                        throw(ArgumentError(
                            "gbc_rhs[$i] contains an AffExpr with a variable $(var) that is not a copied master variable. " *
                            "gbc_rhs should only reference copied master variables passed to customize(). "
                        ))
                    end
                end
            else
                throw(ArgumentError(
                            "gbc_rhs[$i] contains a non-affine expression $rhs. " *
                            "we currently only support affine `gbc_rhs`. "
                        ))
            end
        end

        # Validate gbc_sense: must be a vector of GBCBoundType
        for (i, sense) in enumerate(gbc_sense)
            if !(sense isa GBCBoundType)
                throw(ArgumentError("gbc_sense[$i] is not a GBCBoundType (got $(sense)). Each entry must be among `UpperBound`, `LowerBound`, and `Fixed`."))
            end
        end
        
        return Vector{VariableRef}(gbc_lhs), Vector{Union{VariableRef, AffExpr}}(gbc_rhs), Vector{GBCBoundType}(gbc_sense)
    end
    
    throw(ArgumentError(
        "Invalid return format of `customize` function. Expected nothing or (gbc_lhs, gbc_rhs, gbc_sense). " *
        "Got $(length(result))."
    ))
end

function generate_cuts(oracle::ClassicalOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    set_time_limit_sec(oracle.model, time_limit)
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    
    # Set GBC bounds based on expression evaluation
    _set_gbc_bounds!(oracle.gbc_lhs, oracle.gbc_rhs, oracle.gbc_sense, x_value)
    
    optimize!(oracle.model)
    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end
    
    status = dual_status(oracle.model)
    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(oracle.model)

        a_x = dual.(oracle.fixed_x_constraints)
        
        # Accumulate GBC dual values
        _accumulate_gbc_duals!(a_x, oracle.gbc_lhs, oracle.gbc_rhs, oracle.gbc_sense)
        
        a_t = [-1.0] 
        a_0 = sub_obj_val - a_x'*x_value 
        if sub_obj_val >= t_value[1] * (1 + oracle.param.rtol) + oracle.param.atol / tol_normalize
            return false, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        else
            return true, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        end

    elseif status == INFEASIBILITY_CERTIFICATE || termination_status(oracle.model) == INFEASIBLE
        # Use L1 normalization for feasibility cuts
        return _generate_l1_feasibility_cut(oracle, x_value, t_value, tol_normalize, time_limit)
    else
        throw(UnexpectedModelStatusException("ClassicalOracle: $(status). This is likely a numerical issue. Please try using other oracles, such as unified oracle or pareto oracle."))
    end
end

"""
    _set_gbc_bounds!(gbc_lhs, gbc_rhs, gbc_sense, x_value)

Set the bounds for GBC LHS variables by evaluating RHS expressions with given x values.
Optimized to handle VariableRef directly without AffExpr overhead.
"""
function _set_gbc_bounds!(gbc_lhs::Vector{VariableRef}, 
                          gbc_rhs::Vector{Union{VariableRef, AffExpr}},
                          gbc_sense::Vector{GBCBoundType},
                          x_value::Vector{Float64})
    for i in 1:length(gbc_lhs)
        rhs = gbc_rhs[i]
        # Evaluate RHS expression
        if rhs isa VariableRef
            # Fast path for simple variable bounds
            bound_value = x_value[rhs.index.value]
        else
            # General path for affine expressions
            bound_value = value(v -> x_value[v.index.value], rhs)
        end
        
        # Set bound based on sense
        if gbc_sense[i] == UpperBound
            set_upper_bound(gbc_lhs[i], bound_value)
        elseif gbc_sense[i] == LowerBound
            set_lower_bound(gbc_lhs[i], bound_value)
        else  # Fixed
            fix(gbc_lhs[i], bound_value; force=true)
        end
    end
end

"""
    _restore_gbc_base_bounds!(gbc_lhs, gbc_base_lower, gbc_base_upper)

Restore GBC LHS variable bounds to their base (non-GBC) bounds.
This is used before solving the L1-normalized feasibility problem so that
GBC bounds do not make the normalized model infeasible.
"""
function _restore_gbc_base_bounds!(gbc_lhs::Vector{VariableRef},
                                  gbc_base_lower::Vector{Float64},
                                  gbc_base_upper::Vector{Float64})
    for i in 1:length(gbc_lhs)
        if isfinite(gbc_base_lower[i])
            set_lower_bound(gbc_lhs[i], gbc_base_lower[i])
        elseif has_lower_bound(gbc_lhs[i])
            delete_lower_bound(gbc_lhs[i])
        end

        if isfinite(gbc_base_upper[i])
            set_upper_bound(gbc_lhs[i], gbc_base_upper[i])
        elseif has_upper_bound(gbc_lhs[i])
            delete_upper_bound(gbc_lhs[i])
        end
    end
end

"""
    _accumulate_gbc_duals!(a_x, gbc_lhs, gbc_rhs, gbc_sense)

Accumulate dual values from GBC constraints into the cut coefficients.
Optimized to handle VariableRef directly without AffExpr overhead.
"""
function _accumulate_gbc_duals!(a_x::Vector{Float64},
                                gbc_lhs::Vector{VariableRef},
                                gbc_rhs::Vector{Union{VariableRef, AffExpr}},
                                gbc_sense::Vector{GBCBoundType})
    for i in 1:length(gbc_lhs)
        # Get dual value based on bound sense
        if gbc_sense[i] == UpperBound
            dual_val = dual(UpperBoundRef(gbc_lhs[i]))
        elseif gbc_sense[i] == LowerBound
            dual_val = dual(LowerBoundRef(gbc_lhs[i]))
        else  # Fixed
            dual_val = dual(FixRef(gbc_lhs[i]))
        end
        
        # Accumulate to corresponding x positions
        rhs = gbc_rhs[i]
        if rhs isa VariableRef
            # Fast path for simple variable bounds
            a_x[rhs.index.value] += dual_val
        else
            # General path for affine expressions
            for (var, coef) in rhs.terms
                a_x[var.index.value] += coef * dual_val
            end
        end
    end
end

# ============================================================================
# L1 Normalization Functions for Feasibility Cuts
# ============================================================================

"""
    _enter_l1_mode!(oracle::ClassicalOracle)

Enter L1 normalization mode by:
1. Activating z0 coefficients in >= and <= constraints
2. Activating mu+/mu- coefficients in == constraints  
3. Setting objective to minimize z0

Per L1.md formulation:
- For >= constraints: add +z0 (relaxes the constraint)
- For <= constraints: add -z0 (relaxes the constraint)
- For == constraints: activate mu+ (coeff +1.0) and mu- (coeff -1.0)
"""
function _enter_l1_mode!(oracle::ClassicalOracle)
    state = oracle.l1_state
    # Activate z0 in >= constraints: coefficient +1.0
    for con in state.l1_geq_constraints
        set_normalized_coefficient(con, state.z0, 1.0)
    end
    
    # Activate z0 in <= constraints: coefficient -1.0
    for con in state.l1_leq_constraints
        set_normalized_coefficient(con, state.z0, -1.0)
    end
    
    # Activate mu+/mu- in == constraints
    for (i, con) in enumerate(state.l1_eq_constraints)
        set_normalized_coefficient(con, state.mu_plus[i], 1.0)
        set_normalized_coefficient(con, state.mu_minus[i], -1.0)
    end
    
    # Set objective to minimize z0
    @objective(oracle.model, Min, state.z0)
end

"""
    _exit_l1_mode!(oracle::ClassicalOracle)

Exit L1 normalization mode by:
1. Deactivating z0 coefficients (set to 0)
2. Deactivating mu+/mu- coefficients (set to 0)
3. Restoring original objective
"""
function _exit_l1_mode!(oracle::ClassicalOracle)
    state = oracle.l1_state
    # Deactivate z0 in >= constraints
    for con in state.l1_geq_constraints
        set_normalized_coefficient(con, state.z0, 0.0)
    end
    
    # Deactivate z0 in <= constraints
    for con in state.l1_leq_constraints
        set_normalized_coefficient(con, state.z0, 0.0)
    end
    
    # Deactivate mu+/mu- in == constraints
    for (i, con) in enumerate(state.l1_eq_constraints)
        set_normalized_coefficient(con, state.mu_plus[i], 0.0)
        set_normalized_coefficient(con, state.mu_minus[i], 0.0)
    end
    
    # Restore original objective
    set_objective_function(oracle.model, state.original_objective)
    set_objective_sense(oracle.model, state.original_objective_sense)
end

"""
    _generate_l1_feasibility_cut(oracle, x_value, t_value, tol_normalize, time_limit)

Generate an L1-normalized feasibility cut when the standard subproblem is infeasible.

Algorithm:
1. Enter L1 mode (activate z0/mu coefficients, set min z0 objective)
2. Solve the normalized problem
3. Extract feasibility cut from dual multipliers
4. Exit L1 mode (restore original model)

Note: GBC bounds are NOT set in L1 mode (excluded from normalization per paper).
Note: GBC duals are NOT accumulated in feasibility cuts.
"""
function _generate_l1_feasibility_cut(oracle::ClassicalOracle, 
                                       x_value::Vector{Float64}, 
                                       t_value::Vector{Float64},
                                       tol_normalize::Float64,
                                       time_limit::Float64)
    state = oracle.l1_state
    # Enter L1 mode
    _enter_l1_mode!(oracle)
    
    # Note: GBC bounds are NOT set in L1 mode per paper.
    # The normalized model excludes GBC from z0 penalty.
    _restore_gbc_base_bounds!(oracle.gbc_lhs, state.gbc_base_lower, state.gbc_base_upper)
    
    # Solve the L1 normalized problem
    set_time_limit_sec(oracle.model, time_limit)
    optimize!(oracle.model)
    
    # Check for time limit
    if termination_status(oracle.model) == TIME_LIMIT
        _exit_l1_mode!(oracle)
        throw(TimeLimitException("Time limit reached during L1 normalized feasibility cut generation"))
    end
    
    term_status = termination_status(oracle.model)
    
    if term_status == OPTIMAL
        z0_val = objective_value(oracle.model)
        
        @debug "L1 normalized feasibility cut: z0* = $z0_val"
        
        # Extract duals from fixing constraints
        a_x = dual.(oracle.fixed_x_constraints)
        
        # Note: GBC duals are NOT accumulated for L1 feasibility cuts
        
        # For feasibility cuts: a_t = 0 (no epigraph term)
        a_t = [0.0]
        
        # By strong duality: z0* = π(r - Tx*), so a_0 = z0* - a_x'*x*
        a_0 = z0_val - a_x' * x_value
        
        # Exit L1 mode (restore)
        _exit_l1_mode!(oracle)
        
        # Cut is violated if z0* > 0
        if z0_val >= oracle.param.zero_tol / tol_normalize
            return false, [Hyperplane(a_x, a_t, a_0)], [Inf]
        else
            return true, [Hyperplane(a_x, a_t, a_0)], [Inf]
        end
    else
        # This shouldn't happen - L1 normalized problem should always be feasible
        # because z0 can grow arbitrarily large
        _exit_l1_mode!(oracle)
        throw(UnexpectedModelStatusException(
            "ClassicalOracle L1 mode: Normalized model returned $(term_status). " *
            "This is unexpected as the normalized problem should always be feasible."
        ))
    end
end
