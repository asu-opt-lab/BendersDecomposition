export auto_decompose, auto_decompose_unified, default_master_solver_params, default_oracle_solver_params

using JuMP
using MathOptInterface
const MOI = MathOptInterface

"""
    default_master_solver_params() -> Dict{String, Any}

Return default solver parameters for master problem.
"""
function default_master_solver_params()
    return Dict(
        "solver" => "CPLEX",
        "CPX_PARAM_EPINT" => 1e-9,
        "CPX_PARAM_EPRHS" => 1e-9,
        "CPX_PARAM_EPGAP" => 1e-9,
        "CPXPARAM_Threads" => 4
    )
end

"""
    default_oracle_solver_params() -> Dict{String, Any}

Return default solver parameters for oracle problem.
"""
function default_oracle_solver_params()
    return Dict(
        "solver" => "CPLEX",
        "CPX_PARAM_EPRHS" => 1e-9,
        "CPX_PARAM_NUMERICALEMPHASIS" => 1,
        "CPX_PARAM_EPOPT" => 1e-9
    )
end

function classify_variables_for_benders(model::Model, decision_vars::Vector{Symbol})
    all_vars = all_variables(model)
    decision_vars_set = Set(decision_vars)

    master_vars = VariableRef[]
    oracle_vars = VariableRef[]

    for var in all_vars
        var_name = name(var)
        var_symbol = Symbol(var_name)

        # Check if explicitly specified as decision variable
        if var_symbol in decision_vars_set
            push!(master_vars, var)
        # Check variable type - binary and integer go to master
        elseif is_binary(var) || is_integer(var)
            push!(master_vars, var)
        # All other variables (continuous) go to oracle
        else
            push!(oracle_vars, var)
        end
    end

    return master_vars, oracle_vars
end

"""
    partition_constraints_for_benders(model::Model, master_vars::Vector{VariableRef}, oracle_vars::Vector{VariableRef}) -> Tuple{Vector{ConstraintRef}, Vector{ConstraintRef}, Vector{ConstraintRef}}

Partition constraints based on which variables they involve.
"""
function partition_constraints_for_benders(model::Model, master_vars::Vector{VariableRef}, oracle_vars::Vector{VariableRef})
    master_constraints = ConstraintRef[]
    oracle_constraints = ConstraintRef[]
    coupling_constraints = ConstraintRef[]

    # Check for unsupported variable-in-set constraints
    unsupported_constraint_types = String[]
    for (F, S) in list_of_constraint_types(model)
        if F == Vector{VariableRef} || F == Vector{AffExpr}
            set_name = string(S)
            if occursin("SecondOrderCone", set_name)
                push!(unsupported_constraint_types, "second-order cone constraints (SOCP)")
            elseif occursin("RotatedSecondOrderCone", set_name)
                push!(unsupported_constraint_types, "rotated second-order cone constraints")
            elseif occursin("PositiveSemidefinite", set_name)
                push!(unsupported_constraint_types, "semidefinite constraints (SDP/PSD)")
            elseif occursin("SOS1", set_name)
                push!(unsupported_constraint_types, "special ordered set type 1 constraints (SOS1)")
            elseif occursin("SOS2", set_name)
                push!(unsupported_constraint_types, "special ordered set type 2 constraints (SOS2)")
            elseif occursin("ExponentialCone", set_name)
                push!(unsupported_constraint_types, "exponential cone constraints")
            elseif occursin("PowerCone", set_name)
                push!(unsupported_constraint_types, "power cone constraints")
            else
                # Generic vector constraint
                push!(unsupported_constraint_types, "vector/conic constraints ($set_name)")
            end
        elseif F == QuadExpr
            push!(unsupported_constraint_types, "quadratic constraints")
        elseif F == NonlinearExpr
            push!(unsupported_constraint_types, "nonlinear constraints")
        end
    end

    if !isempty(unsupported_constraint_types)
        unique_types = unique(unsupported_constraint_types)
        throw(ArgumentError("auto_decompose does not support the following constraint types found in the model: " *
                          join(unique_types, ", ") * ". " *
                          "Only linear constraints (affine expressions) are supported."))
    end

    # Get all constraints
    model_constraints = all_constraints(model, include_variable_in_set_constraints=false)

    master_vars_set = Set(master_vars)
    oracle_vars_set = Set(oracle_vars)

    for constraint in model_constraints
        constraint_vars = Set(constraint_variables(constraint))

        has_master_vars = !isempty(intersect(constraint_vars, master_vars_set))
        has_oracle_vars = !isempty(intersect(constraint_vars, oracle_vars_set))

        if has_master_vars && has_oracle_vars
            # Coupling constraint
            push!(coupling_constraints, constraint)
        elseif has_master_vars
            # Pure master constraint
            push!(master_constraints, constraint)
        elseif has_oracle_vars
            # Pure oracle constraint
            push!(oracle_constraints, constraint)
        end
    end

    return master_constraints, oracle_constraints, coupling_constraints
end

"""
    decompose_objective_for_benders(model::Model, master_vars::Vector{VariableRef}, oracle_vars::Vector{VariableRef}) -> Tuple{AffExpr, AffExpr}

Decompose the objective function into master and oracle parts.
"""
function decompose_objective_for_benders(model::Model, master_vars::Vector{VariableRef}, oracle_vars::Vector{VariableRef})
    obj = objective_function(model)

    master_vars_set = Set(master_vars)
    oracle_vars_set = Set(oracle_vars)

    master_terms = AffExpr(0.0)
    oracle_terms = AffExpr(0.0)

    # Handle constant term
    if isa(obj, AffExpr)
        master_terms.constant = obj.constant

        # Process linear terms
        for (var, coeff) in obj.terms
            if var in master_vars_set
                add_to_expression!(master_terms, coeff, var)
            elseif var in oracle_vars_set
                add_to_expression!(oracle_terms, coeff, var)
            end
        end
    elseif isa(obj, VariableRef)
        # Single variable objective
        if obj in master_vars_set
            add_to_expression!(master_terms, 1.0, obj)
        elseif obj in oracle_vars_set
            add_to_expression!(oracle_terms, 1.0, obj)
        end
    else
        throw(ArgumentError("Unsupported objective function type in auto_decompose: $(typeof(obj)). " *
                          "Only linear objectives (AffExpr) and single variable objectives (VariableRef) are supported."))
    end

    return master_terms, oracle_terms
end

"""
    set_variable_bounds!(new_var::VariableRef, original_var::VariableRef)

Set bounds on new variable based on original variable bounds.
"""
function set_variable_bounds!(new_var::VariableRef, original_var::VariableRef)
    if is_fixed(original_var)
        fix(new_var, fix_value(original_var); force=true)
        return
    end
    if has_lower_bound(original_var)
        set_lower_bound(new_var, lower_bound(original_var))
    end
    if has_upper_bound(original_var)
        set_upper_bound(new_var, upper_bound(original_var))
    end
end

"""
    copy_variable_with_attributes(var::VariableRef, target_model::Model) -> VariableRef

Copy a variable with all its attributes to the target model.
"""
function copy_variable_with_attributes(var::VariableRef, target_model::Model)
    var_name = name(var)

    if is_binary(var)
        new_var = @variable(target_model, binary=true, base_name=var_name)
    elseif is_integer(var)
        new_var = @variable(target_model, integer=true, base_name=var_name)
        set_variable_bounds!(new_var, var)
    else
        new_var = @variable(target_model, base_name=var_name)
        set_variable_bounds!(new_var, var)
    end

    return new_var
end

"""
    copy_variables_to_model(vars::Vector{VariableRef}, target_model::Model) -> Tuple{Vector{VariableRef}, Dict{VariableRef, VariableRef}}

Copy multiple variables to target model with all attributes and return new variables and mapping.
"""
function copy_variables_to_model(vars::Vector{VariableRef}, target_model::Model)
    var_map = Dict{VariableRef, VariableRef}()
    new_vars = VariableRef[]

    for var in vars
        new_var = copy_variable_with_attributes(var, target_model)
        var_map[var] = new_var
        push!(new_vars, new_var)
    end

    return new_vars, var_map
end

"""
    copy_constraints_to_model(constraints::Vector{ConstraintRef}, target_model::Model, var_map::Dict{VariableRef, VariableRef})

Copy multiple constraints to target model with variable substitution using provided mapping.
"""
function copy_constraints_to_model(constraints::Vector{ConstraintRef}, target_model::Model, var_map::Dict{VariableRef, VariableRef})
    for constraint in constraints
        copy_constraint_with_substitution(constraint, target_model, var_map)
    end
end

"""
    substitute_variables_in_expression(expr, var_map::Dict{VariableRef, VariableRef})

Substitute variables in an expression using the provided mapping.
"""
function substitute_variables_in_expression(expr, var_map::Dict{VariableRef, VariableRef})
    if isa(expr, AffExpr)
        new_expr = AffExpr(expr.constant)
        for (var, coeff) in expr.terms
            add_to_expression!(new_expr, coeff, var_map[var])
        end
        return new_expr
    elseif isa(expr, VariableRef)
        return var_map[expr]
    else
        throw(ArgumentError("Unsupported expression type in auto_decompose: $(typeof(expr)). " *
                          "Only linear expressions (AffExpr) and variable references (VariableRef) are supported."))
    end
end

"""
    copy_constraint_with_substitution(constraint::ConstraintRef, target_model::Model, var_map::Dict{VariableRef, VariableRef}) -> ConstraintRef

Copy a constraint to the target model with variable substitution.
"""
function copy_constraint_with_substitution(constraint::ConstraintRef, target_model::Model, var_map::Dict{VariableRef, VariableRef})
    constraint_obj = constraint_object(constraint)
    substituted_func = substitute_variables_in_expression(constraint_obj.func, var_map)
    substituted_set = constraint_obj.set
    return add_constraint(target_model, ScalarConstraint(substituted_func, substituted_set))
end

"""
    constraint_variables(constraint::ConstraintRef) -> Vector{VariableRef}

Get variables involved in a constraint.
"""
function constraint_variables(constraint::ConstraintRef)
    constraint_obj = constraint_object(constraint)
    func = constraint_obj.func

    if isa(func, AffExpr)
        return collect(keys(func.terms))
    elseif isa(func, VariableRef)
        return [func]
    else
        throw(ArgumentError("Unsupported constraint function type in auto_decompose: $(typeof(func)). " *
                          "Only linear constraints are supported. " *
                          "The following constraint types are NOT supported: " *
                          "quadratic constraints, nonlinear constraints, " *
                          "second-order cone constraints (SOCP), " *
                          "semidefinite constraints (SDP/PSD), " *
                          "special ordered set constraints (SOS1/SOS2), " *
                          "and other conic/vector constraints."))
    end
end

"""
    extract_variable_coefficients(objective::AffExpr, vars::Vector{VariableRef}) -> Vector{Float64}

Extract coefficients for specified variables from the objective function.
"""
function extract_variable_coefficients(objective::AffExpr, vars::Vector{VariableRef})
    coeffs = zeros(length(vars))
    for (i, var) in enumerate(vars)
        if haskey(objective.terms, var)
            coeffs[i] = objective.terms[var]
        end
    end
    return coeffs
end

"""
    build_data_from_decomposition(master_vars::Vector{VariableRef}, master_objective::AffExpr) -> Data

Build Data structure from decomposition information.
"""
function build_data_from_decomposition(master_vars::Vector{VariableRef}, master_objective::AffExpr)
    dim_x = length(master_vars)
    dim_t = 1  # Single auxiliary variable

    # Extract coefficients
    c_x = extract_variable_coefficients(master_objective, master_vars)
    c_t = [1.0]  # Coefficient for auxiliary variable

    # Create simplified problem data
    problem_data = AutoDecomposedData()

    return Data(dim_x, dim_t, problem_data, c_x, c_t)
end

"""
    build_complete_master(master_vars::Vector{VariableRef}, master_constraints::Vector{ConstraintRef},
                         master_objective::AffExpr, original_model::Model, master_solver_param::Dict{String,Any}) -> Master

Build complete Master problem.
"""
function build_complete_master(master_vars::Vector{VariableRef}, master_constraints::Vector{ConstraintRef},
                              master_objective::AffExpr, original_model::Model, master_solver_param::Dict{String,Any})
    master_model = Model()
    assign_attributes!(master_model, master_solver_param)
    set_silent(master_model)

    master_var_refs, var_map = copy_variables_to_model(master_vars, master_model)
    copy_constraints_to_model(master_constraints, master_model, var_map)

    t = @variable(master_model, [1:1], base_name="t", lower_bound=-1e6)

    if !iszero(master_objective)
        substituted_obj = substitute_variables_in_expression(master_objective, var_map)
        @objective(master_model, objective_sense(original_model), substituted_obj + t[1])
    else
        @objective(master_model, objective_sense(original_model), t[1])
    end

    return Master(master_model, master_var_refs, t)
end

"""
    build_complete_classical_oracle(oracle_vars::Vector{VariableRef}, master_vars::Vector{VariableRef},
                                   oracle_constraints::Vector{ConstraintRef}, coupling_constraints::Vector{ConstraintRef},
                                   oracle_objective::AffExpr, original_model::Model, oracle_solver_param::Dict{String,Any}, oracle_param::BasicOracleParam) -> ClassicalOracle

Build complete ClassicalOracle.
"""
function build_complete_classical_oracle(oracle_vars::Vector{VariableRef}, master_vars::Vector{VariableRef},
                                        oracle_constraints::Vector{ConstraintRef}, coupling_constraints::Vector{ConstraintRef},
                                        oracle_objective::AffExpr, original_model::Model, oracle_solver_param::Dict{String,Any}, oracle_param::BasicOracleParam)
    oracle_model = Model()
    assign_attributes!(oracle_model, oracle_solver_param)
    set_silent(oracle_model)

    _, var_map = copy_variables_to_model(oracle_vars, oracle_model)

    # Do not move fixing constriants of master vars to sub
    for var in master_vars
        new_var = @variable(oracle_model, base_name = name(var))
        if is_fixed(var)
            fix(new_var, fix_value(var); force = true)
        end
        var_map[var] = new_var
    end

    copy_constraints_to_model(oracle_constraints, oracle_model, var_map)
    copy_constraints_to_model(coupling_constraints, oracle_model, var_map)

    fixing_constraints = ConstraintRef[]
    for (i, var) in enumerate(master_vars)
        oracle_var = var_map[var]
        fixing_constraint = @constraint(oracle_model, oracle_var == 0)
        push!(fixing_constraints, fixing_constraint)
    end

    if !iszero(oracle_objective)
        substituted_obj = substitute_variables_in_expression(oracle_objective, var_map)
        @objective(oracle_model, objective_sense(original_model), substituted_obj)
    end

    return ClassicalOracle(oracle_param, oracle_model, fixing_constraints)
end

function build_complete_unified_oracle(oracle_vars::Vector{VariableRef}, master_vars::Vector{VariableRef},
                                       oracle_constraints::Vector{ConstraintRef}, coupling_constraints::Vector{ConstraintRef},
                                       oracle_objective::AffExpr, original_model::Model,
                                       oracle_solver_param::Dict{String,Any},
                                       oracle_param::BasicOracleParam)

    oracle_model = Model()
    assign_attributes!(oracle_model, oracle_solver_param)
    set_silent(oracle_model)

    # 1) Copy oracle vars
    _, var_map = copy_variables_to_model(oracle_vars, oracle_model)

    # 2) Create oracle-side copies of master vars (x for reformulation)
    x_oracle = VariableRef[]
    for var in master_vars
        v = @variable(oracle_model, base_name = name(var))  # Do not move fixing constriants of master vars to sub
        if is_fixed(var)
            fix(v, fix_value(var); force = true)
        end
        var_map[var] = v
        push!(x_oracle, v)
    end

    # 3) Copy constraints
    copy_constraints_to_model(oracle_constraints, oracle_model, var_map)
    copy_constraints_to_model(coupling_constraints, oracle_model, var_map)

    # 4) Set oracle objective before reformulation
    if !iszero(oracle_objective)
        substituted_obj = substitute_variables_in_expression(oracle_objective, var_map)
        @objective(oracle_model, objective_sense(original_model), substituted_obj)
    else
        @objective(oracle_model, objective_sense(original_model), 0.0)
    end

    # 5) Reformulate in-place (w0 is hard-coded to 1.0 inside)
    model_reformulation!(oracle_model; x = x_oracle)

    # 6) Return unified oracle
    return UnifiedOracle(oracle_param, oracle_model)
end

"""
    auto_decompose(model::Model; decision_vars::Vector{Symbol} = Symbol[], master_solver_param::Dict{String,Any} = default_master_solver_params(), oracle_solver_param::Dict{String,Any} = default_oracle_solver_params(), oracle_param::BasicOracleParam = BasicOracleParam()) -> Tuple{Data, Master, ClassicalOracle}

Automatically decompose a JuMP model into Data, Master, and ClassicalOracle components.

# Limitations
- Only supports linear (affine) objective functions. Quadratic and nonlinear objectives are not supported.
- Only supports linear constraints. Quadratic, nonlinear, and vector constraints are not supported.
- The model must be decomposable with clear first-stage (binary/integer) and second-stage (continuous) variable structure.

# Arguments
- `model::Model`: The JuMP model to decompose
- `decision_vars::Vector{Symbol}`: List of variable names to be treated as master variables (default: empty, uses automatic classification)
- `master_solver_param::Dict{String,Any}`: Solver parameters for the master problem (default: default_master_solver_params())
- `oracle_solver_param::Dict{String,Any}`: Solver parameters for the oracle problem (default: default_oracle_solver_params())
- `oracle_param::BasicOracleParam`: Oracle-specific parameters (default: BasicOracleParam())

# Returns
- `Tuple{Data, Master, ClassicalOracle}`: The decomposed problem components

# Examples
```julia
# Basic usage
data, master, oracle = auto_decompose(model)

# With custom solver parameters
data, master, oracle = auto_decompose(model;
    master_solver_param = Dict("solver" => "CPLEX", "CPXPARAM_Threads" => 8),
    oracle_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPOPT" => 1e-12)
)

# With custom oracle parameters
data, master, oracle = auto_decompose(model;
    oracle_param = BasicOracleParam(rtol = 1e-6, atol = 1e-6)
)
```
"""
function auto_decompose(model::Model;
                       decision_vars::Vector{Symbol} = Symbol[],
                       master_solver_param::Dict{String,Any} = default_master_solver_params(),
                       oracle_solver_param::Dict{String,Any} = default_oracle_solver_params(),
                       oracle_param::BasicOracleParam = BasicOracleParam())

    # Check for empty model
    if num_variables(model) == 0
        throw(ArgumentError("Cannot decompose empty model with no variables in auto_decompose."))
    end

    # Check for unsupported objective types
    obj = objective_function(model)
    if isa(obj, QuadExpr) || isa(obj, NonlinearExpr)
        throw(ArgumentError("Quadratic and nonlinear objectives are not supported in auto_decompose."))
    end

    # Step 1: Classify variables
    master_vars, oracle_vars = classify_variables_for_benders(model, decision_vars)

    # Step 2: Partition constraints
    master_constraints, oracle_constraints, coupling_constraints = partition_constraints_for_benders(model, master_vars, oracle_vars)

    # Step 3: Decompose objective function
    master_objective, oracle_objective = decompose_objective_for_benders(model, master_vars, oracle_vars)

    # Step 4: Build Data structure
    data = build_data_from_decomposition(master_vars, master_objective)

    # Step 5: Build Master problem
    master = build_complete_master(master_vars, master_constraints, master_objective, model, master_solver_param)

    # Step 6: Build ClassicalOracle
    oracle = build_complete_classical_oracle(oracle_vars, master_vars, oracle_constraints, coupling_constraints,
                                            oracle_objective, model, oracle_solver_param, oracle_param)

    return data, master, oracle
end

function auto_decompose_unified(model::Model;
                                decision_vars::Vector{Symbol} = Symbol[],
                                master_solver_param::Dict{String,Any} = default_master_solver_params(),
                                oracle_solver_param::Dict{String,Any} = default_oracle_solver_params(),
                                oracle_param::BasicOracleParam = BasicOracleParam())

    if num_variables(model) == 0
        throw(ArgumentError("Cannot decompose empty model with no variables."))
    end

    obj = objective_function(model)
    if isa(obj, QuadExpr) || isa(obj, NonlinearExpr)
        throw(ArgumentError("Quadratic and nonlinear objectives are not supported."))
    end

    master_vars, oracle_vars = classify_variables_for_benders(model, decision_vars)

    master_constraints, oracle_constraints, coupling_constraints =
        partition_constraints_for_benders(model, master_vars, oracle_vars)

    master_objective, oracle_objective =
        decompose_objective_for_benders(model, master_vars, oracle_vars)

    data = build_data_from_decomposition(master_vars, master_objective)

    master = build_complete_master(master_vars, master_constraints, master_objective, model, master_solver_param)

    oracle = build_complete_unified_oracle(oracle_vars, master_vars,
                                           oracle_constraints, coupling_constraints,
                                           oracle_objective, model,
                                           oracle_solver_param, oracle_param)

    return data, master, oracle
end

function model_reformulation!(model::Model; x)
    w0 = 1.0

    # Weight in dual problem
    @variable(model, z)

    # Add epigraph constraint
    @constraint(model, epigraph, w0 * z .- objective_function(model) .>= 0)

    # Change objective function
    @objective(model, Min, z)

    # Linking constraints
    @constraint(model, fix_x_lb, z .+ x .>= 0)
    @constraint(model, fix_x_ub, z .- x .>= 0)
end