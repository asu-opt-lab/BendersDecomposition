export ParetoOracle, ParetoOracleParam

"""
    ParetoOracleParam <: AbstractOracleParam

Parameter structure for ParetoOracle implementing Magnanti-Wong Pareto-optimal cuts.

# Fields
- `rtol::Float64`: Relative tolerance for cut violation detection (default: 1e-9)
- `atol::Float64`: Absolute tolerance for cut violation detection (default: 0.0)
- `zero_tol::Float64`: Threshold below which a value is considered zero (default: 1e-9)
- `core_point::Vector{Float64}`: Core point x_0 for Magnanti-Wong problem (REQUIRED).
"""
struct ParetoOracleParam <: AbstractOracleParam
    rtol::Float64
    atol::Float64
    zero_tol::Float64
    core_point::Vector{Float64}

    function ParetoOracleParam(core_point::Vector{Float64}; rtol = 1e-9, atol = 0.0, zero_tol = 1e-9)
        isempty(core_point) && throw(ArgumentError("core_point must be provided and non-empty"))
        new(rtol, atol, zero_tol, core_point)
    end
end

"""
    ParetoOracle <: AbstractTypicalOracle

A Pareto oracle for Benders decomposition that generates Pareto-optimal cuts using 
the Magnanti-Wong technique.

The Pareto oracle maintains two models:
1. `model`: Standard subproblem (same as ClassicalOracle)
2. `pareto_model`: Magnanti-Wong primal problem for Pareto-optimal cuts

## Magnanti-Wong Primal Problem
```math
\\begin{align*}
\\min \\quad & d^{\\top}y + \\xi^*z \\\\
\\text{s.t.} \\quad & Ax + By + bz \\geq b \\quad (\\pi_0) \\\\
& x + x^*z = x_0 \\quad (\\pi_1) \\\\
& y \\geq 0
\\end{align*}
```

Where:
- `ξ*` is the optimal objective value from the standard subproblem
- `x*` is the current master solution
- `x_0` is the core point (parameter)
- The cut coefficients come from π_1* (duals of pareto_fixing_constraints)

# Fields
- `param::ParetoOracleParam`: Oracle parameters including core_point
- `model::Model`: Standard subproblem model
- `fixed_x_constraints::Vector{ConstraintRef}`: Fixing constraints in standard model
- `pareto_model::Model`: Magnanti-Wong primal problem model
- `pareto_variable::VariableRef`: The σ (z) variable in pareto_model
- `pareto_fixing_constraints::Vector{ConstraintRef}`: Fixing constraints in pareto_model (x + x*σ = x_0)

# Constructor
```julia
ParetoOracle(data::AbstractData, master::Master, param::ParetoOracleParam;
             customize = customize_sub_model!,
             scen_idx::Int = 0)
```
"""
mutable struct ParetoOracle <: AbstractTypicalOracle
    
    param::ParetoOracleParam

    # Standard subproblem model
    model::Model
    fixed_x_constraints::Vector{ConstraintRef}

    # Magnanti-Wong pareto model
    pareto_model::Model
    pareto_variable::VariableRef
    pareto_fixing_constraints::Vector{ConstraintRef}

    function ParetoOracle(data::AbstractData, master::Master, param::ParetoOracleParam;
                         customize = customize_sub_model!,
                         scen_idx::Int = 0)
    
        @debug "Building Pareto oracle"
        model = Model()

        # Copy the master's coupling variables into the submodel (with identical axes and symbols)
        x_copy = copy_variables!(model, master.x_tuple)

        # Collect all copied master variables
        x = var_from_tuple(x_copy)

        # Validate core_point dimension
        dim_x = length(x)
        if length(param.core_point) != dim_x
            throw(DimensionMismatch("core_point has length $(length(param.core_point)) but expected $dim_x"))
        end

        # Build the submodel using user-defined customization
        # NOTE: Do NOT add fixing constraints yet
        customize(model, data, scen_idx; x_copy...)

        # Validate that all constraints are supported types (LP)
        _validate_constraint_types(model)

        # ---------------------------------------------------------
        # Build Pareto Model (Re-run construction instead of copying)
        # ---------------------------------------------------------
        pareto_model = Model()
        
        # Copy master variables for pareto model
        pareto_x_copy = copy_variables!(pareto_model, master.x_tuple)
        pareto_x = var_from_tuple(pareto_x_copy)
        
        # Build pareto model structure
        customize(pareto_model, data, scen_idx; pareto_x_copy...)
        
        # Apply Magnanti-Wong transformations (add σ, etc.)
        pareto_variable, pareto_fixing_constraints = _apply_pareto_transformations!(pareto_model, pareto_x)

        # ---------------------------------------------------------
        # Finalize Standard Model
        # ---------------------------------------------------------
        # NOW add fixing constraints to standard model
        @constraint(model, fix_x, x .== 0)

        new(param, model, fix_x,
            pareto_model, pareto_variable, pareto_fixing_constraints)
    end

    ParetoOracle() = new()
end

# Add constructor for SeparableOracle compatibility
function ParetoOracle(data::AbstractData, master::Master; 
                      customize = customize_sub_model!,
                      scen_idx::Int = 0,
                      param::ParetoOracleParam)
    return ParetoOracle(data, master, param; customize = customize, scen_idx = scen_idx)
end


"""
    _validate_constraint_types(model::Model)

Validate that all constraints in the model are supported types for ParetoOracle.
Only GreaterThan (>=), LessThan (<=), and EqualTo (==) constraints are supported.
Throws UnsupportedModelException if unsupported constraint types are found.
"""
function _validate_constraint_types(model::Model)
    # 1. Check for Integer/Binary variables (Mixed-Integer terms)
    for (F, S) in list_of_constraint_types(model)
        if S <: Union{MOI.Integer, MOI.ZeroOne, MOI.Semicontinuous, MOI.Semiinteger}
            throw(UnsupportedModelException(
                "Unsupported constraint type: $S. " *
                "ParetoOracle requires a continuous Linear Programming (LP) subproblem. " *
                "Integer or binary variables are not allowed."
            ))
        end
    end

    # 2. Check Objective Function (must be linear)
    obj_type = objective_function_type(model)
    if obj_type != Nothing && !(obj_type <: Union{VariableRef, AffExpr, Real})
        throw(UnsupportedModelException(
            "Unsupported objective function type: $obj_type. " *
            "ParetoOracle only supports Linear Programming (LP) with linear objectives."
        ))
    end

    # 3. Check Structural Constraints (Linearity and Supported Sets)
    for con in all_constraints(model, include_variable_in_set_constraints=false)
        con_obj = constraint_object(con)
        set = con_obj.set
        func = con_obj.func
        
        # Check constraint set type (EqualTo, GreaterThan, LessThan)
        if !(set isa MOI.GreaterThan || set isa MOI.LessThan || set isa MOI.EqualTo)
            throw(UnsupportedModelException(
                "Unsupported constraint set type: $(typeof(set)). " *
                "ParetoOracle only supports GreaterThan (>=), LessThan (<=), and EqualTo (==) constraints. " *
                "If you have interval or other constraint types, please reformulate them."
            ))
        end

        # Check constraint function type (must be linear: AffExpr or VariableRef)
        if !(func isa AffExpr || func isa VariableRef)
            throw(UnsupportedModelException(
                "Unsupported constraint function type: $(typeof(func)). " *
                "ParetoOracle only supports linear constraints (LP). " *
                "Quadratic or other non-linear constraints are not supported."
            ))
        end
    end
end

"""
    _apply_pareto_transformations!(pareto_model::Model, x_vars::Vector{VariableRef})

Transform a standard subproblem model into a Magnanti-Wong primal problem.

1. Adds σ variable
2. Adds b*σ term to all problem constraints
3. Creates fixing constraints: x + x*σ = x_0

Returns: (pareto_variable, pareto_fixing_constraints)
"""
function _apply_pareto_transformations!(pareto_model::Model, x_vars::Vector{VariableRef})
    
    # Step 1: Add σ variable (the z in Magnanti-Wong formulation)
    # We add a large but finite bound to σ to improve numerical stability
    σ = @variable(pareto_model, σ)
    
    # Step 2: Add b*σ term to all problem constraints
    # Note: Constraint types are validated by _validate_constraint_types on the base model
    for con in all_constraints(pareto_model, include_variable_in_set_constraints=false)
        con_obj = constraint_object(con)
        set = con_obj.set
        rhs = normalized_rhs(con)
        
        if set isa MOI.GreaterThan
            # Ax + By >= b  -->  Ax + By + b*σ >= b
            set_normalized_coefficient(con, σ, rhs)
        elseif set isa MOI.LessThan
            # Ax + By <= b  -->  Ax + By + b*σ <= b
            set_normalized_coefficient(con, σ, rhs)
        else  # MOI.EqualTo
            # Ax + By = b  -->  Ax + By + b*σ = b
            set_normalized_coefficient(con, σ, rhs)
        end
    end
    
    # Step 3: Create fixing constraints for pareto_model: x + x*σ = x_0
    pareto_fixing_constraints = ConstraintRef[]
    for x_var in x_vars
        # x_var is already the variable in pareto_model, no mapping needed!
        con = @constraint(pareto_model, x_var + 0.0 * σ == 0)
        push!(pareto_fixing_constraints, con)
    end
    
    # Step 4: Add σ to objective (coefficient will be set to ξ* in generate_cuts)
    original_obj = objective_function(pareto_model)
    @objective(pareto_model, Min, original_obj + 0.0 * σ)
    
    return σ, pareto_fixing_constraints
end

"""
    generate_cuts(oracle::ParetoOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; 
                  tol_normalize = 1.0, time_limit = 3600)

Generate Pareto-optimal Benders cuts using the Magnanti-Wong technique.

## Algorithm (from ref/cuts.jl):
1. Set x = x* in standard model and solve to get ξ*
2. If feasible:
   - Set objective coefficient of σ to ξ*
   - Set σ coefficients in fixing constraints to x*
   - Set RHS of fixing constraints to core_point x_0
   - Solve pareto_model
   - Get cut coefficients from duals of pareto_fixing_constraints
3. If infeasibility certificate:
   - Use classical feasibility cut from standard model duals

## Returns
- `is_in_L::Bool`: Whether the point is in the feasible region L
- `hyperplanes::Vector{Hyperplane}`: Generated cuts
- `sub_obj_vals::Vector{Float64}`: Subproblem objective values
"""
function generate_cuts(oracle::ParetoOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; 
                       tol_normalize = 1.0, time_limit = 3600)
    
    # Set time limits
    set_time_limit_sec(oracle.model, time_limit)
    set_time_limit_sec(oracle.pareto_model, time_limit)
    
    # Step 1: Set x = x* in standard model 
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    
    # Step 2: Solve standard model to get ξ*
    optimize!(oracle.model)
    
    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during Pareto cut generation (standard model)"))
    end
    
    status = dual_status(oracle.model)
    
    if status == FEASIBLE_POINT
        # Get optimal objective value ξ* from standard model
        sub_obj_val = objective_value(oracle.model)
        
        # Step 3: Set up pareto_model for Magnanti-Wong problem
        # Set objective coefficient of σ to ξ*
        set_objective_coefficient(oracle.pareto_model, oracle.pareto_variable, sub_obj_val)
        
        # Set σ coefficient in fixing constraints to x*
        # Constraint: x + x*·σ = x_0
        # We filter small coefficients to improve numerical scaling
        for i in 1:length(x_value)
            coef = abs(x_value[i]) > oracle.param.zero_tol ? x_value[i] : 0.0
            set_normalized_coefficient(oracle.pareto_fixing_constraints[i], oracle.pareto_variable, coef)
        end

        # Set RHS to core_point x_0 
        set_normalized_rhs.(oracle.pareto_fixing_constraints, oracle.param.core_point)
        
        # Step 4: Solve pareto_model
        optimize!(oracle.pareto_model)
        
        if termination_status(oracle.pareto_model) == TIME_LIMIT
            throw(TimeLimitException("Time limit reached during Pareto cut generation (pareto model)"))
        end
        
        pareto_status = dual_status(oracle.pareto_model)
        
        if pareto_status == FEASIBLE_POINT || pareto_status == MOI.NEARLY_FEASIBLE_POINT
            # Get cut coefficients from pareto_fixing_constraints duals
            a_x = dual.(oracle.pareto_fixing_constraints)
            
            # Cut: t >= ξ* - π_1*'x* + π_1*'x
            a_t = [-1.0]
            a_0 = sub_obj_val - dot(a_x, x_value)
            
            # Check if cut is violated
            if sub_obj_val >= t_value[1] * (1 + oracle.param.rtol) + oracle.param.atol / tol_normalize
                return false, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
            else
                return true, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
            end
        else
            @warn "ParetoOracle: Unexpected dual status $(pareto_status). This is likely a numerical issue. Falling back to typical cut."
            a_x = dual.(oracle.fixed_x_constraints)
            a_t = [-1.0]
            a_0 = sub_obj_val - dot(a_x, x_value)
            
            if sub_obj_val >= t_value[1] * (1 + oracle.param.rtol) + oracle.param.atol / tol_normalize
                return false, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
            else
                return true, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
            end
        end
        
    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(oracle.model)
            dual_sub_obj_val = dual_objective_value(oracle.model)
            @debug "dual_sub_obj_val = $dual_sub_obj_val"
            
            a_x = dual.(oracle.fixed_x_constraints)
            a_t = [0.0]
            a_0 = dual_sub_obj_val - dot(a_x, x_value)
            
            if dual_sub_obj_val >= oracle.param.zero_tol / tol_normalize
                return false, [Hyperplane(a_x, a_t, a_0)], [Inf]
            else
                return true, [Hyperplane(a_x, a_t, a_0)], [Inf]
            end
        end
    else
        throw(UnexpectedModelStatusException("ParetoOracle: Unexpected dual status $(status). This is likely a numerical issue."))
    end
end
