"""
    L1NormOracle - L1 Normalization for CGLP (Cut Generating Linear Program)

Based on Section 4.2 of "Implementing Automatic Benders Decomposition in a Modern MIP Solver"

The textbook Benders CGLP gives little control over which feasibility cut is returned.
To select "good" feasibility cuts, we add L1 normalization that truncates the dual cone.

## Two-Stage Approach:
1. First solve standard CGLP: if feasible → derive optimality cut
2. If infeasible → add z0 penalty variable and solve normalized problem
3. The z0 acts as L1-norm normalization in dual space, favoring cuts with sparse support

## L1 Normalization (equations 17-19 from paper):
When standard subproblem is infeasible:
    min  z0
    s.t. Tx* + Qy + z0 >= r  (add z0 to all constraints involving master vars)
         y, z0 >= 0

By strong duality: z0* = π(r - Tx*), hence the cut is violated by x*.

## Important Notes:
- GBC (Generalized Bound Constraints) are NOT included in normalization
- If d=0 (no objective), never remove normalization
"""

# Import from BendersX (provides all BendersBase + BendersLibrary exports)
using BendersX
using JuMP
using LinearAlgebra
using SparseArrays
const MOI = JuMP.MOI  # Access MOI through JuMP

# CRITICAL: Import generate_cuts to extend it (not shadow it)
import BendersX: generate_cuts

export L1NormOracle, L1NormOracleParam

"""
    L1NormOracleParam <: AbstractOracleParam

Parameter structure for L1NormOracle.

# Fields
- `rtol::Float64`: Relative tolerance for cut violation detection (default: 1e-9)
- `atol::Float64`: Absolute tolerance for cut violation detection (default: 0.0)
- `zero_tol::Float64`: Threshold below which a value is considered zero (default: 1e-9)
"""
struct L1NormOracleParam <: AbstractOracleParam
    rtol::Float64
    atol::Float64
    zero_tol::Float64

    function L1NormOracleParam(; rtol = 1e-9, atol = 0.0, zero_tol = 1e-9)
        new(rtol, atol, zero_tol)
    end
end

"""
    L1NormOracle <: AbstractTypicalOracle

An oracle for Benders decomposition that uses L1 normalization for feasibility cuts.

This oracle maintains two models:
1. `model`: Standard subproblem (same as ClassicalOracle)
2. `norm_model`: Normalized model with z0 penalty variable for feasibility cuts

## Algorithm:
1. Solve standard subproblem with x fixed to x*
2. If feasible → return optimality cut (same as ClassicalOracle)
3. If infeasible → solve normalized model with z0 penalty → return L1-normalized feasibility cut

# Fields
- `param::L1NormOracleParam`: Oracle parameters
- `model::Model`: Standard subproblem model
- `fixed_x_constraints::Vector{ConstraintRef}`: Fixing constraints in standard model
- `gbc_lhs::Vector{VariableRef}`: GBC left-hand-side variables
- `gbc_rhs::Vector{Union{VariableRef, AffExpr}}`: GBC right-hand-side expressions
- `gbc_sense::Vector{GBCBoundType}`: GBC bound types
- `norm_model::Model`: Normalized model with z0 variable
- `norm_fixed_x_constraints::Vector{ConstraintRef}`: Fixing constraints in normalized model
- `norm_z0::VariableRef`: The z0 penalty variable
- `has_objective::Bool`: Whether the original subproblem has a non-zero objective
"""
mutable struct L1NormOracle <: AbstractTypicalOracle
    
    param::L1NormOracleParam

    # Standard subproblem model (same as ClassicalOracle)
    model::Model
    fixed_x_constraints::Vector{ConstraintRef}

    # GBC (Generalized Bound Constraints) for standard model
    gbc_lhs::Vector{VariableRef}
    gbc_rhs::Vector{Union{VariableRef, AffExpr}}
    gbc_sense::Vector{GBCBoundType}

    # Normalized model for feasibility cuts
    norm_model::Model
    norm_fixed_x_constraints::Vector{ConstraintRef}
    norm_z0::VariableRef
    
    # GBC for normalized model (NOT included in z0 normalization, but still need bounds/duals)
    norm_gbc_lhs::Vector{VariableRef}
    norm_gbc_rhs::Vector{Union{VariableRef, AffExpr}}
    norm_gbc_sense::Vector{GBCBoundType}
    
    # Track if original objective is zero (d=0 case)
    has_objective::Bool

    function L1NormOracle(data::AbstractData, master::Master; 
                          customize = customize_sub_model!,
                          scen_idx::Int = 0, 
                          param::L1NormOracleParam = L1NormOracleParam())
    
        @debug "Building L1 normalization oracle"
        
        # ========================================
        # Build Standard Model (same as Classical)
        # ========================================
        model = Model()

        # Copy the master's coupling variables into the submodel
        x_copy = copy_variables!(model, master.x_tuple)
        x = var_from_tuple(x_copy)
        
        # Add linking constraint
        @constraint(model, fix_x, x .== 0)

        # Build the submodel using user-defined customization
        result = customize(model, data, scen_idx; x_copy...)
        
        # Validate constraint types
        _validate_constraint_types(model)
        
        # Parse GBC result
        gbc_lhs, gbc_rhs, gbc_sense = _parse_gbc_result(result, x)
        
        # Check if objective is non-zero
        obj_func = objective_function(model)
        has_objective = !_is_zero_objective(obj_func)

        # ========================================
        # Build Normalized Model with z0
        # ========================================
        norm_model = Model()
        
        # Copy master variables for normalized model
        norm_x_copy = copy_variables!(norm_model, master.x_tuple)
        norm_x = var_from_tuple(norm_x_copy)
        
        # Add linking constraint for normalized model
        @constraint(norm_model, norm_fix_x, norm_x .== 0)
        
        # Build normalized submodel structure
        norm_result = customize(norm_model, data, scen_idx; norm_x_copy...)
        
        # Parse GBC result for normalized model
        norm_gbc_lhs, norm_gbc_rhs, norm_gbc_sense = _parse_gbc_result(norm_result, norm_x)
        
        # Get set of fixing constraint indices to exclude from normalization
        fix_x_indices = Set(c.index for c in norm_fix_x)
        
        # Add z0 penalty variable and apply L1 normalization
        norm_z0 = _apply_l1_normalization!(norm_model, norm_x, fix_x_indices)

        new(param, model, fix_x, gbc_lhs, gbc_rhs, gbc_sense,
            norm_model, norm_fix_x, norm_z0, 
            norm_gbc_lhs, norm_gbc_rhs, norm_gbc_sense,
            has_objective)
    end

    L1NormOracle() = new()
end

"""
    _is_zero_objective(obj_func) -> Bool

Check if the objective function is effectively zero (d=0 case).
"""
function _is_zero_objective(obj_func)
    if obj_func isa Real
        return abs(obj_func) < 1e-10
    elseif obj_func isa VariableRef
        return false  # Variable means non-zero objective
    elseif obj_func isa AffExpr
        # Check if all coefficients are zero
        if abs(obj_func.constant) > 1e-10
            return false
        end
        for (_, coef) in obj_func.terms
            if abs(coef) > 1e-10
                return false
            end
        end
        return true
    end
    return false
end

"""
    _apply_l1_normalization!(norm_model::Model, fixed_x_vars::Vector{VariableRef}, excluded_constraint_indices::Set)

Apply L1 normalization to the model by adding z0 penalty variable.

L1 Normalization (equations 17-19 from paper):
    min  z0
    s.t. Tx* + Qy + z0 >= r  (for >= constraints)
         Tx* + Qy - z0 <= r  (for <= constraints)
         y, z0 >= 0

The z0 is added to constraints that involve the master (x) variables,
acting as normalization in the dual space.

Note: GBC constraints are NOT modified (they don't contain x variables in the constraint body).
Note: Constraints in `excluded_constraint_indices` are NOT modified (e.g., fixing constraints).

Returns: z0 variable reference
"""
function _apply_l1_normalization!(norm_model::Model, fixed_x_vars::Vector{VariableRef}, excluded_constraint_indices::Set)
    # Add z0 penalty variable (z0 >= 0)
    z0 = @variable(norm_model, z0 >= 0)
    
    # Build set of master variable indices for fast lookup
    x_indices = Set(v.index for v in fixed_x_vars)
    
    # Collect constraints to process first (avoid modifying during iteration)
    constraints_to_process = ConstraintRef[]
    for con in all_constraints(norm_model, include_variable_in_set_constraints=false)
        # Skip excluded constraints (e.g., fixing constraints)
        if con.index in excluded_constraint_indices
            continue
        end
        if _constraint_involves_x(con, x_indices)
            push!(constraints_to_process, con)
        end
    end
    
    # Add z0 to constraints that involve master variables
    for con in constraints_to_process
        _add_z0_to_constraint!(norm_model, con, z0)
    end
    
    # Set objective to minimize z0
    @objective(norm_model, Min, z0)
    
    return z0
end

"""
    _constraint_involves_x(con::ConstraintRef, x_indices::Set) -> Bool

Check if a constraint involves any of the master (x) variables.
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
    _add_z0_to_constraint!(model::Model, con::ConstraintRef, z0::VariableRef)

Add z0 to a constraint based on its sense:
- For >= constraints: add +z0 to LHS (makes constraint easier to satisfy)
- For <= constraints: add -z0 to LHS (makes constraint easier to satisfy)
- For == constraints: split into >= and <=, each with z0

This ensures z0* represents the minimum violation needed.
"""
function _add_z0_to_constraint!(model::Model, con::ConstraintRef, z0::VariableRef)
    con_obj = constraint_object(con)
    func = con_obj.func
    set = con_obj.set
    
    if set isa MOI.GreaterThan
        # f(x,y) >= r  -->  f(x,y) + z0 >= r
        set_normalized_coefficient(con, z0, 1.0)
    elseif set isa MOI.LessThan
        # f(x,y) <= r  -->  f(x,y) - z0 <= r
        set_normalized_coefficient(con, z0, -1.0)
    elseif set isa MOI.EqualTo
        # == constraint: split into >= and <= with z0
        # f(x,y) = r  -->  f(x,y) + z0 >= r AND f(x,y) - z0 <= r
        rhs = set.value
        original_name = name(con)
        
        # Create >= constraint: f(x) + z0 >= rhs
        lb_con = @constraint(model, func >= rhs)
        set_normalized_coefficient(lb_con, z0, 1.0)
        if !isempty(original_name)
            set_name(lb_con, original_name * "_lb")
        end
        
        # Create <= constraint: f(x) - z0 <= rhs
        ub_con = @constraint(model, func <= rhs)
        set_normalized_coefficient(ub_con, z0, -1.0)
        if !isempty(original_name)
            set_name(ub_con, original_name * "_ub")
        end
        
        # Delete original equality constraint
        delete(model, con)
    end
end

"""
    generate_cuts(oracle::L1NormOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; 
                  tol_normalize = 1.0, time_limit = 3600)

Generate Benders cuts using L1 normalization for feasibility cuts.

## Two-Stage Algorithm:
1. Set x = x* in standard model and solve
2. If feasible (dual status = FEASIBLE_POINT):
   - Extract duals and return optimality cut (same as ClassicalOracle)
3. If infeasible (dual status = INFEASIBILITY_CERTIFICATE):
   - Set x = x* in normalized model
   - Solve min z0 subject to normalized constraints
   - Extract duals and return L1-normalized feasibility cut

## Returns
- `is_in_L::Bool`: Whether the point is in the feasible region L
- `hyperplanes::Vector{Hyperplane}`: Generated cuts
- `sub_obj_vals::Vector{Float64}`: Subproblem objective values
"""
function generate_cuts(oracle::L1NormOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; 
                                   tol_normalize = 1.0, time_limit = 3600)
    
    set_time_limit_sec(oracle.model, time_limit)
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    
    # Set GBC bounds based on expression evaluation
    _set_gbc_bounds!(oracle.gbc_lhs, oracle.gbc_rhs, oracle.gbc_sense, x_value)
    
    optimize!(oracle.model)
    
    term_status = termination_status(oracle.model)
    
    if term_status == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during L1Norm cut generation"))
    end
    
    status = dual_status(oracle.model)
    
    if status == FEASIBLE_POINT
        # ========================================
        # Stage 1: Optimality Cut (same as Classical)
        # ========================================
        sub_obj_val = objective_value(oracle.model)

        a_x = dual.(oracle.fixed_x_constraints)
        
        # Accumulate GBC dual values
        _accumulate_gbc_duals!(a_x, oracle.gbc_lhs, oracle.gbc_rhs, oracle.gbc_sense)
        
        a_t = [-1.0] 
        a_0 = sub_obj_val - a_x' * x_value 
        
        if sub_obj_val >= t_value[1] * (1 + oracle.param.rtol) + oracle.param.atol / tol_normalize
            return false, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        else
            return true, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        end

    elseif status == INFEASIBILITY_CERTIFICATE || term_status == INFEASIBLE
        # ========================================
        # Stage 2: L1 Normalized Feasibility Cut
        # ========================================
        # For L1NormOracle, we don't need the Farkas certificate from the standard model.
        # When infeasible, we switch to the normalized model which is always feasible
        # (z0 can grow arbitrarily large to satisfy all constraints).
        return _generate_l1_normalized_feasibility_cut(oracle, x_value, t_value, tol_normalize, time_limit)
    else
        throw(UnexpectedModelStatusException("L1NormOracle: term_status=$(term_status), dual_status=$(status). This is likely a numerical issue."))
    end
end

"""
    _generate_l1_normalized_feasibility_cut(oracle, x_value, t_value, tol_normalize, time_limit)

Generate an L1-normalized feasibility cut when the standard subproblem is infeasible.

This solves the normalized problem:
    min  z0
    s.t. Tx* + Qy + z0 >= r
         y, z0 >= 0

And extracts the cut from the dual multipliers.

Note: GBC bounds are NOT set in the normalized model. GBC are excluded from normalization
per the paper, meaning they don't participate in the z0 penalty. This allows the normalized
model to remain feasible even when GBC would make it infeasible. However, GBC duals are NOT
accumulated since the bounds aren't set.
"""
function _generate_l1_normalized_feasibility_cut(oracle::L1NormOracle, 
                                                  x_value::Vector{Float64}, 
                                                  t_value::Vector{Float64},
                                                  tol_normalize::Float64,
                                                  time_limit::Float64)
    
    set_time_limit_sec(oracle.norm_model, time_limit)
    set_normalized_rhs.(oracle.norm_fixed_x_constraints, x_value)
    
    # Note: GBC bounds are NOT set in normalized model.
    # GBC constraints are excluded from normalization per paper - they don't get z0 penalty.
    # Setting GBC bounds here would make the problem infeasible when standard subproblem is infeasible.
    # The feasibility cut from the normalized model doesn't include GBC contributions.
    
    optimize!(oracle.norm_model)
    
    if termination_status(oracle.norm_model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during L1 normalized feasibility cut generation"))
    end
    
    norm_status = termination_status(oracle.norm_model)
    
    if norm_status == OPTIMAL
        z0_val = objective_value(oracle.norm_model)
        
        @debug "L1 normalized feasibility cut: z0* = $z0_val"
        
        # Extract duals from fixing constraints
        a_x = dual.(oracle.norm_fixed_x_constraints)
        
        # Note: GBC duals are NOT accumulated since GBC bounds are not set in normalized model
        
        # For feasibility cuts: a_t = 0 (no epigraph term)
        a_t = [0.0]
        
        # By strong duality: z0* = π(r - Tx*), so a_0 = z0* - a_x'*x*
        a_0 = z0_val - a_x' * x_value
        
        # Cut is violated if z0* > 0
        if z0_val >= oracle.param.zero_tol / tol_normalize
            return false, [Hyperplane(a_x, a_t, a_0)], [Inf]
        else
            return true, [Hyperplane(a_x, a_t, a_0)], [Inf]
        end
    else
        # This shouldn't happen - normalized problem should always be feasible
        # because z0 can grow arbitrarily large
        throw(UnexpectedModelStatusException(
            "L1NormOracle: Normalized model returned $(norm_status). " *
            "This is unexpected as the normalized problem should always be feasible."
        ))
    end
end

# ============================================================================
# Helper functions from ClassicalOracle (need to be accessible)
# ============================================================================

"""
    _set_gbc_bounds!(gbc_lhs, gbc_rhs, gbc_sense, x_value)

Set the bounds for GBC LHS variables by evaluating RHS expressions with given x values.
"""
function _set_gbc_bounds!(gbc_lhs::Vector{VariableRef}, 
                          gbc_rhs::Vector{Union{VariableRef, AffExpr}},
                          gbc_sense::Vector{GBCBoundType},
                          x_value::Vector{Float64})
    for i in 1:length(gbc_lhs)
        rhs = gbc_rhs[i]
        # Evaluate RHS expression
        if rhs isa VariableRef
            bound_value = x_value[rhs.index.value]
        else
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
    _accumulate_gbc_duals!(a_x, gbc_lhs, gbc_rhs, gbc_sense)

Accumulate dual values from GBC constraints into the cut coefficients.
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
            a_x[rhs.index.value] += dual_val
        else
            for (var, coef) in rhs.terms
                a_x[var.index.value] += coef * dual_val
            end
        end
    end
end

"""
    _validate_constraint_types(model::Model)

Validate that all constraints in the model are supported types for typical oracles.
"""
function _validate_constraint_types(model::Model)
    # Check for Integer/Binary variables
    for (F, S) in list_of_constraint_types(model)
        if S <: Union{MOI.Integer, MOI.ZeroOne, MOI.Semicontinuous, MOI.Semiinteger}
            throw(UnsupportedModelException(
                "Unsupported constraint type: $S. " *
                "Typical oracles require a continuous LP subproblem."
            ))
        end
    end

    # Check Objective Function
    obj_type = objective_function_type(model)
    if obj_type != Nothing && !(obj_type <: Union{VariableRef, AffExpr, Real})
        throw(UnsupportedModelException(
            "Unsupported objective function type: $obj_type. " *
            "Typical oracles only support linear objectives."
        ))
    end

    # Check Structural Constraints
    for con in all_constraints(model, include_variable_in_set_constraints=false)
        con_obj = constraint_object(con)
        set = con_obj.set
        func = con_obj.func
        
        if !(set isa MOI.GreaterThan || set isa MOI.LessThan || set isa MOI.EqualTo)
            throw(UnsupportedModelException(
                "Unsupported constraint set type: $(typeof(set)). " *
                "Only >=, <=, and == constraints are supported."
            ))
        end

        if !(func isa AffExpr || func isa VariableRef)
            throw(UnsupportedModelException(
                "Unsupported constraint function type: $(typeof(func)). " *
                "Only linear constraints are supported."
            ))
        end
    end
end

"""
    _parse_gbc_result(result, x_vars) -> (gbc_lhs, gbc_rhs, gbc_sense)

Parse the value returned by a user `customize` function for GBC.
"""
function _parse_gbc_result(result, x_vars::Vector{VariableRef})
    if result === nothing
        return VariableRef[], Union{VariableRef, AffExpr}[], GBCBoundType[]
    end
    
    if length(result) == 3
        gbc_lhs, gbc_rhs, gbc_sense = result
        
        n = length(gbc_lhs)
        if length(gbc_rhs) != n || length(gbc_sense) != n
            throw(DimensionMismatch("All GBC vectors must have the same length."))
        end
        
        if isempty(gbc_lhs)
            return VariableRef[], Union{VariableRef, AffExpr}[], GBCBoundType[]
        end
        
        # Validation (simplified from ClassicalOracle)
        x_indices = Set(v.index for v in x_vars)
        
        for (i, lhs_var) in enumerate(gbc_lhs)
            if !(lhs_var isa VariableRef)
                throw(ArgumentError("gbc_lhs[$i] is not a VariableRef."))
            end
            if lhs_var.index in x_indices
                throw(ArgumentError("gbc_lhs[$i] contains a copied master variable."))
            end
        end
        
        for (i, rhs) in enumerate(gbc_rhs)
            if rhs isa VariableRef
                if !(rhs.index in x_indices)
                    throw(ArgumentError("gbc_rhs[$i] contains a non-master variable."))
                end
            elseif rhs isa AffExpr
                for (var, _) in rhs.terms
                    if !(var.index in x_indices)
                        throw(ArgumentError("gbc_rhs[$i] AffExpr contains a non-master variable."))
                    end
                end
            else
                throw(ArgumentError("gbc_rhs[$i] is not a VariableRef or AffExpr."))
            end
        end
        
        for (i, sense) in enumerate(gbc_sense)
            if !(sense isa GBCBoundType)
                throw(ArgumentError("gbc_sense[$i] is not a GBCBoundType."))
            end
        end
        
        return Vector{VariableRef}(gbc_lhs), Vector{Union{VariableRef, AffExpr}}(gbc_rhs), Vector{GBCBoundType}(gbc_sense)
    end
    
    throw(ArgumentError("Invalid return format of `customize` function."))
end
