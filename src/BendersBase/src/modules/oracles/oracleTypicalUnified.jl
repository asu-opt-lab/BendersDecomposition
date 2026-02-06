export UnifiedOracle, UnifiedOracleParam

"""
    UnifiedOracleParam <: AbstractOracleParam

Parameter structure for UnifiedOracle.

# Fields
- `rtol::Float64`: Relative tolerance for cut violation detection (default: 1e-9)
- `atol::Float64`: Absolute tolerance for cut violation detection (default: 0.0)
- `zero_tol::Float64`: Threshold below which a value is considered zero (default: 1e-9)
- `w0::Float64`: Weight for the epigraph constraint in the unified formulation (default: 1.0).
  Controls the relative importance of objective violation vs constraint violation.
"""
struct UnifiedOracleParam <: AbstractOracleParam
    rtol::Float64
    atol::Float64
    zero_tol::Float64
    w0::Float64

    function UnifiedOracleParam(; rtol = 1e-9, atol = 0.0, zero_tol = 1e-9, w0 = 1.0)
        w0 > 0 || throw(ArgumentError("UnifiedOracleParam: w0 must be positive, got $w0"))
        new(rtol, atol, zero_tol, w0)
    end
end

"""
    UnifiedOracle <: AbstractTypicalOracle

An oracle that generates unified Benders cuts via a sigma-relaxation reformulation.

UnifiedOracle reformulates the classical subproblem by:
1. Adding a nonnegative slack variable 
2. Relaxing constraints that involve master variables based on constraint sense
3. Replacing each fixing equality with lower- and upper-bound fixing constraints
4. Converting the original objective bound into a constraint with weight `w0`

This approach generates stronger cuts and can handle both feasible and infeasible 
subproblems in a unified manner.

# Fields
- `param::UnifiedOracleParam`: Oracle parameters (includes w0 weight)
- `model::Model`: JuMP model in unified form
- `fixing_lb_constraints::Vector{ConstraintRef}`: Lower-bound fixing constraints in `model`
- `fixing_ub_constraints::Vector{ConstraintRef}`: Upper-bound fixing constraints in `model`
- `objective_constraint::ConstraintRef`: Objective-bound constraint in `model`

# Constructor
```julia
UnifiedOracle(data::AbstractData, master::Master; 
              customize = customize_sub_model!,
              scen_idx::Int = 0, 
              param::UnifiedOracleParam = UnifiedOracleParam())
```

# Example with custom w0
```julia
param = UnifiedOracleParam(w0 = 2.0)  # Higher weight for objective violation
oracle = UnifiedOracle(data, master; param = param)
```
"""
mutable struct UnifiedOracle <: AbstractTypicalOracle
    
    param::UnifiedOracleParam

    model::Model
    
    # Unified constraint structure
    fixing_lb_constraints::Vector{ConstraintRef}
    fixing_ub_constraints::Vector{ConstraintRef}
    objective_constraint::ConstraintRef

    function UnifiedOracle(data::AbstractData, master::Master; 
                          customize = customize_sub_model!,
                          scen_idx::Int = 0, 
                          param::UnifiedOracleParam = UnifiedOracleParam())
    
        @debug "Building unified oracle"
        model = Model()

        # Copy the master's coupling variables into the submodel (with identical axes and symbols)
        x_copy = copy_variables!(model, master.x_tuple)

        # Collect all copied master variables
        x = var_from_tuple(x_copy)

        # Build the submodel using user-defined customization, passing the copied variables
        customize(model, data, scen_idx; x_copy...)

        # Validate that all constraints are supported types (LP)
        _validate_constraint_types(model)

        # Create oracle instance
        oracle = new()
        oracle.param = param
        oracle.model = model
        oracle.fixing_lb_constraints = ConstraintRef[]
        oracle.fixing_ub_constraints = ConstraintRef[]
        
        # Transform to unified form (x is passed for decision variable detection)
        _apply_unified_transformations!(oracle, x)
        
        return oracle
    end

    UnifiedOracle() = new()
end

"""
    _apply_unified_transformations!(oracle::UnifiedOracle, fixed_x_vars::Vector{VariableRef})

Transform the classical subproblem model into the unified form with a sigma relaxation variable.
"""
function _apply_unified_transformations!(oracle::UnifiedOracle, fixed_x_vars::Vector{VariableRef})
    model = oracle.model
    w0 = oracle.param.w0
    
    # Step 1: Get original objective function (must be done before constraint changes)
    original_objective = objective_function(model)
    
    # Step 2: Add nonnegative sigma variable
    σ = @variable(model, σ >= 0)
    
    # Step 3: Rewrite problem constraints with sigma BEFORE creating new constraints
    # This eliminates the need to track which constraints to skip
    _rewrite_problem_constraints_with_sigma!(model, σ, fixed_x_vars)
    
    # Step 4: Create unified fixing constraints 
    # Original: x == value  -->  x + σ >= value AND x - σ <= value
    for x_var in fixed_x_vars
        lb_con = @constraint(model, x_var + σ >= 0)
        push!(oracle.fixing_lb_constraints, lb_con)
        
        ub_con = @constraint(model, x_var - σ <= 0)
        push!(oracle.fixing_ub_constraints, ub_con)
    end
    
    # Step 5: Convert original objective to constraint with w0 weight
    # For minimization: original_obj <= t  -->  -original_obj + w0*σ >= -t
    oracle.objective_constraint = @constraint(model, -original_objective + w0 * σ >= 0)
    
    # Step 6: Set new objective to minimize σ
    @objective(model, Min, σ)
end

"""
    _rewrite_problem_constraints_with_sigma!(model::Model, σ::VariableRef, fixed_x_vars::Vector{VariableRef})

Rewrite problem constraints to include sigma (`σ`) based on constraint sense.
Only constraints that contain decision variables (fixed_x_vars) are relaxed.
"""
function _rewrite_problem_constraints_with_sigma!(model::Model, σ::VariableRef, fixed_x_vars::Vector{VariableRef})
    # Build set of decision variables for fast lookup
    decision_vars_set = Set{VariableRef}(fixed_x_vars)
    
    # Process each constraint (excluding variable bounds)
    for con in all_constraints(model, include_variable_in_set_constraints=false)
        if _contains_decision_vars(con, decision_vars_set)
            _rewrite_single_constraint!(model, con, σ)
        end
    end
end

"""
    _contains_decision_vars(con::ConstraintRef, decision_vars_set::Set{VariableRef}) -> Bool

Check if a constraint contains any of the decision variables.
Returns true if the constraint's function references any variable in the decision set.

Note: Constraint function type validation is performed by _validate_constraint_types() 
before this function is called, so we only handle VariableRef and AffExpr.
"""
function _contains_decision_vars(con::ConstraintRef, decision_vars_set::Set{VariableRef})
    func = constraint_object(con).func
    
    if func isa VariableRef
        return func in decision_vars_set
    else  # AffExpr (validated by _validate_constraint_types)
        return any(var in decision_vars_set for (var, _) in func.terms)
    end
end

"""
    _insert_suffix(name::String, suffix::String) -> String

Insert suffix before the index brackets in a constraint name.
Examples:
- `demand[1]` + `_lb` -> `demand_lb[1]`
- `flow[i,j]` + `_ub` -> `flow_ub[i,j]`
- `simple` + `_lb` -> `simple_lb`
"""
function _insert_suffix(name::String, suffix::String)
    bracket_pos = findfirst('[', name)
    if bracket_pos === nothing
        return name * suffix
    else
        return name[1:bracket_pos-1] * suffix * name[bracket_pos:end]
    end
end

"""
    _rewrite_single_constraint!(model::Model, con::ConstraintRef, σ::VariableRef)

Relax the constraints by adding slack variable.
"""
function _rewrite_single_constraint!(model::Model, con::ConstraintRef, σ::VariableRef)
    
    # Get constraint object
    con_obj = constraint_object(con)
    func = con_obj.func
    set = con_obj.set
    
    if set isa MOI.GreaterThan
        # >= constraint: add +σ coefficient
        # Original: f(x) >= rhs  -->  f(x) + σ >= rhs
        set_normalized_coefficient(con, σ, 1.0)
        
    elseif set isa MOI.LessThan
        # <= constraint: add -σ coefficient
        # Original: f(x) <= rhs  -->  f(x) - σ <= rhs
        set_normalized_coefficient(con, σ, -1.0)
        
    else  # MOI.EqualTo (validated by _validate_constraint_types)
        # == constraint: split into >= and <= with σ
        rhs = set.value
        original_name = name(con)
        
        # Create >= constraint: f(x) + σ >= rhs
        lb_con = @constraint(model, func >= rhs)
        set_normalized_coefficient(lb_con, σ, 1.0)
        if !isempty(original_name)
            set_name(lb_con, _insert_suffix(original_name, "_lb"))
        end
        
        # Create <= constraint: f(x) - σ <= rhs
        ub_con = @constraint(model, func <= rhs)
        set_normalized_coefficient(ub_con, σ, -1.0)
        if !isempty(original_name)
            set_name(ub_con, _insert_suffix(original_name, "_ub"))
        end
        
        # Delete original equality constraint
        delete(model, con)
    end
end

"""
    generate_cuts(oracle::UnifiedOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; 
                  tol_normalize = 1.0, time_limit = 3600)

Generate Benders cuts using the unified oracle.
"""
function generate_cuts(oracle::UnifiedOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; 
                       tol_normalize = 1.0, time_limit = 3600)
    set_time_limit_sec(oracle.model, time_limit)
    
    set_normalized_rhs(oracle.objective_constraint, -t_value[1])
    set_normalized_rhs.(oracle.fixing_lb_constraints, x_value)
    set_normalized_rhs.(oracle.fixing_ub_constraints, x_value)
    
    optimize!(oracle.model)
    
    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during unified cut generation"))
    elseif termination_status(oracle.model) !== OPTIMAL
        throw(UnexpectedModelStatusException("UnifiedOracle: Unexpected termination status $(termination_status(oracle.model))."))
    end
    
    status = dual_status(oracle.model)
    
    if status == FEASIBLE_POINT

        σ_val = objective_value(oracle.model)

        decision_coeffs = dual.(oracle.fixing_lb_constraints) .+ dual.(oracle.fixing_ub_constraints)
        auxiliary_coeffs = dual(oracle.objective_constraint)
        const_term = σ_val + auxiliary_coeffs * t_value[1] - dot(decision_coeffs, x_value)
        
        a_x = decision_coeffs
        a_t = [-auxiliary_coeffs]
        a_0 = const_term
        
        if abs(σ_val) <= oracle.param.zero_tol
            return true, [Hyperplane(a_x, a_t, a_0)], [t_value[1]]
        end
        return false, [Hyperplane(a_x, a_t, a_0)], [Inf]
    else
        throw(UnexpectedModelStatusException("UnifiedOracle: Unexpected dual status $(status)."))
    end
end
