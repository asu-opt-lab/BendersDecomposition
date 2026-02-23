export auto_decompose, default_dcglp_solver_params, default_dcglp_param, default_disjunctive_oracle_param

# Import auto_decompose to extend it
import BendersBase: auto_decompose

# Import internal helper functions from BendersBase
import BendersBase: classify_variables_for_benders, partition_constraints_for_benders,
                     decompose_objective_for_benders, build_data_from_decomposition,
                     build_complete_master, build_complete_classical_oracle, build_complete_unified_oracle,
                     substitute_variables_in_expression, constraint_variables

using JuMP
using MathOptInterface
const MOI = MathOptInterface

"""
    default_dcglp_solver_params() -> Dict{String, Any}

Return default solver parameters for DCGLP model in DisjunctiveOracle.
"""
function default_dcglp_solver_params()
    return Dict(
        "solver" => "CPLEX",
        "CPX_PARAM_EPRHS" => 1e-9,
        "CPX_PARAM_NUMERICALEMPHASIS" => 1,
        "CPX_PARAM_EPOPT" => 1e-9
    )
end

"""
    default_dcglp_param() -> DcglpParam

Return default DCGLP algorithm parameters for DisjunctiveOracle.
"""
function default_dcglp_param()
    return DcglpParam(
        time_limit = 1000.0,
        gap_tolerance = 1e-3,
        halt_limit = 3,
        iter_limit = 250,
        verbose = false
    )
end

"""
    default_disjunctive_oracle_param() -> DisjunctiveOracleParam

Return default DisjunctiveOracle parameters.
"""
function default_disjunctive_oracle_param()
    return DisjunctiveOracleParam(
        norm = LpNorm(Inf),
        split_index_selection_rule = RandomFractional(),
        disjunctive_cut_append_rule = AllDisjunctiveCuts(),
        strengthened = true,
        add_benders_cuts_to_master = 1,
        fraction_of_benders_cuts_to_master = 1.0,
        reuse_dcglp = true,
        lift = false,
        adjust_t_to_fx = false
    )
end

"""
    add_master_constraints_to_dcglp!(dcglp::Model, master_constraints::Vector{ConstraintRef}, master_vars::Vector{VariableRef})

Add master constraints to DCGLP model in disjunctive form.
For each master constraint a'x >= b, creates constraints for i=1:2:
    sum(a[j] * omega_x[i,j]) >= b * omega_0[i]
"""
function add_master_constraints_to_dcglp!(dcglp::Model, master_constraints::Vector{ConstraintRef}, master_vars::Vector{VariableRef})
    # Return if there are no master constraints
    if isempty(master_constraints)
        return
    end

    # Create mapping from variables to indices
    var_to_index = Dict(var => i for (i, var) in enumerate(master_vars))

    omega_x = dcglp[:omega_x]
    omega_0 = dcglp[:omega_0]

    # Iterate through each master constraint
    for constraint in master_constraints
        constraint_obj = constraint_object(constraint)
        func = constraint_obj.func  # AffExpr
        set = constraint_obj.set

        # Extract coefficients and constant term
        # func is in the form: a'x + c
        constant_term = func.constant  # c

        # Build coefficient vector
        coeffs = zeros(length(master_vars))
        for (var, coeff) in func.terms
            if haskey(var_to_index, var)
                idx = var_to_index[var]
                coeffs[idx] = coeff
            else
                throw(ArgumentError("Master constraint contains variable not in master_vars"))
            end
        end

        # Process based on constraint type
        if isa(set, MOI.GreaterThan)
            # a'x + c >= b  =>  a'x >= b - c
            rhs = set.lower - constant_term
            for i in 1:2
                @constraint(dcglp, sum(coeffs[j] * omega_x[i,j] for j in 1:length(master_vars)) >= rhs * omega_0[i])
            end
        elseif isa(set, MOI.LessThan)
            # a'x + c <= b  =>  a'x <= b - c
            rhs = set.upper - constant_term
            for i in 1:2
                @constraint(dcglp, sum(coeffs[j] * omega_x[i,j] for j in 1:length(master_vars)) <= rhs * omega_0[i])
            end
        elseif isa(set, MOI.EqualTo)
            # a'x + c == b  =>  a'x == b - c
            rhs = set.value - constant_term
            for i in 1:2
                @constraint(dcglp, sum(coeffs[j] * omega_x[i,j] for j in 1:length(master_vars)) == rhs * omega_0[i])
            end
        else
            throw(ArgumentError("Unsupported constraint type: $(typeof(set))"))
        end
    end
end

"""
    build_complete_disjunctive_oracle(oracle_vars::Vector{VariableRef}, master_vars::Vector{VariableRef},
                                     master_constraints::Vector{ConstraintRef},
                                     oracle_constraints::Vector{ConstraintRef}, coupling_constraints::Vector{ConstraintRef},
                                     oracle_objective::AffExpr, original_model::Model, data::Data,
                                     typical_oracle_solver_param::Dict{String,Any}, dcglp_solver_param::Dict{String,Any},
                                     oracle_param::DisjunctiveOracleParam, dcglp_param::DcglpParam) -> DisjunctiveOracle

Build complete DisjunctiveOracle with two underlying ClassicalOracles.
"""
function build_complete_disjunctive_oracle(oracle_vars::Vector{VariableRef}, master_vars::Vector{VariableRef},
                                          master_constraints::Vector{ConstraintRef},
                                          oracle_constraints::Vector{ConstraintRef}, coupling_constraints::Vector{ConstraintRef},
                                          oracle_objective::AffExpr, original_model::Model, data::Data,
                                          typical_oracle_solver_param::Dict{String,Any}, dcglp_solver_param::Dict{String,Any},
                                          oracle_param::DisjunctiveOracleParam, dcglp_param::DcglpParam)
    # Build two ClassicalOracles for kappa and nu
    # typical_oracles = [
    #     build_complete_classical_oracle(oracle_vars, master_vars, oracle_constraints, coupling_constraints,
    #                                    oracle_objective, original_model, typical_oracle_solver_param, BasicOracleParam()),
    #     build_complete_classical_oracle(oracle_vars, master_vars, oracle_constraints, coupling_constraints,
    #                                    oracle_objective, original_model, typical_oracle_solver_param, BasicOracleParam())
    # ]

    # Build two UnifiedOracles for kappa and nu
    typical_oracles = [
        build_complete_unified_oracle(oracle_vars, master_vars, oracle_constraints, coupling_constraints,
                                       oracle_objective, original_model, typical_oracle_solver_param, BasicOracleParam()),
        build_complete_unified_oracle(oracle_vars, master_vars, oracle_constraints, coupling_constraints,
                                       oracle_objective, original_model, typical_oracle_solver_param, BasicOracleParam())
    ]

    # Build DisjunctiveOracle
    disjunctive_oracle = DisjunctiveOracle(data, typical_oracles;
                                          solver_param = dcglp_solver_param,
                                          param = dcglp_param)

    # Set oracle-specific parameters
    set_parameter!(disjunctive_oracle, oracle_param)

    # Add master constraints to dcglp in disjunctive form
    add_master_constraints_to_dcglp!(disjunctive_oracle.dcglp, master_constraints, master_vars)

    return disjunctive_oracle
end

"""
    auto_decompose(model::Model, oracle_type::Symbol; decision_vars::Vector{Symbol} = Symbol[], master_solver_param::Dict{String,Any} = default_master_solver_params(), typical_oracle_solver_param::Dict{String,Any} = default_oracle_solver_params(), dcglp_solver_param::Dict{String,Any} = default_dcglp_solver_params(), oracle_param::DisjunctiveOracleParam = default_disjunctive_oracle_param(), dcglp_param::DcglpParam = default_dcglp_param()) -> Tuple{Data, Master, DisjunctiveOracle}

Automatically decompose a JuMP model into Data, Master, and DisjunctiveOracle components.

This is an extension of BendersBase.auto_decompose that supports disjunctive oracle creation.

# Arguments
- `model::Model`: The JuMP model to decompose
- `oracle_type::Symbol`: Type of oracle to create (must be :disjunctive for this method)
- `decision_vars::Vector{Symbol}`: List of variable names to be treated as master variables (default: empty, uses automatic classification)
- `master_solver_param::Dict{String,Any}`: Solver parameters for the master problem (default: default_master_solver_params())
- `typical_oracle_solver_param::Dict{String,Any}`: Solver parameters for typical oracles (default: default_oracle_solver_params())
- `dcglp_solver_param::Dict{String,Any}`: Solver parameters for DCGLP model (default: default_dcglp_solver_params())
- `oracle_param::DisjunctiveOracleParam`: Oracle-specific parameters (default: default_disjunctive_oracle_param())
- `dcglp_param::DcglpParam`: DCGLP algorithm parameters (default: default_dcglp_param())

# Returns
- `Tuple{Data, Master, DisjunctiveOracle}`: The decomposed problem components

# Examples
```julia
# Disjunctive oracle with defaults
data, master, oracle = auto_decompose(model, :disjunctive)

# Disjunctive oracle with custom parameters
data, master, oracle = auto_decompose(model, :disjunctive;
    oracle_param = DisjunctiveOracleParam(strengthened = false),
    dcglp_param = DcglpParam(time_limit = 500.0)
)

# Disjunctive oracle with custom solver parameters
data, master, oracle = auto_decompose(model, :disjunctive;
    master_solver_param = Dict("solver" => "CPLEX", "CPXPARAM_Threads" => 8),
    dcglp_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPOPT" => 1e-12)
)
```
"""
function auto_decompose(model::Model,
                       oracle_type::Symbol;
                       decision_vars::Vector{Symbol} = Symbol[],
                       master_solver_param::Dict{String,Any} = default_master_solver_params(),
                       typical_oracle_solver_param::Dict{String,Any} = default_oracle_solver_params(),
                       dcglp_solver_param::Dict{String,Any} = default_dcglp_solver_params(),
                       oracle_param::DisjunctiveOracleParam = default_disjunctive_oracle_param(),
                       dcglp_param::DcglpParam = default_dcglp_param())
    # Validate oracle_type parameter
    if oracle_type != :disjunctive
        throw(ArgumentError("This method only supports :disjunctive oracle type, got :$oracle_type. For :classical oracle, use BendersBase.auto_decompose without oracle_type parameter."))
    end

    # Check for empty model
    if num_variables(model) == 0
        throw(ArgumentError("Cannot decompose empty model with no variables"))
    end

    # Check for unsupported objective types
    obj = objective_function(model)
    if isa(obj, QuadExpr) || isa(obj, NonlinearExpr)
        throw(ArgumentError("Quadratic and nonlinear objectives are not supported"))
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

    # Step 6: Build DisjunctiveOracle
    oracle = build_complete_disjunctive_oracle(oracle_vars, master_vars, master_constraints,
                                              oracle_constraints, coupling_constraints,
                                              oracle_objective, model, data, typical_oracle_solver_param,
                                              dcglp_solver_param, oracle_param, dcglp_param)

    return data, master, oracle
end
