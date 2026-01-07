export ClassicalOracle, ClassicalOracleParam

"""
    ClassicalOracleParam

Parameter type for [`ClassicalOracle`](@ref).

`ClassicalOracleParam` is an alias of [`BasicOracleParam`](@ref) and controls
numerical tolerances used during cut generation, such as
relative and absolute tolerances for optimality and feasibility checks.

See also: [`ClassicalOracle`](@ref), [`BasicOracleParam`](@ref)
"""
const ClassicalOracleParam = BasicOracleParam

"""
    ClassicalOracle <: AbstractTypicalOracle

Classical Benders cut generator based on dual information from a subproblem.

`ClassicalOracle` implements the standard (classical) Benders decomposition
oracle. Given candidate values of the master variables, it solves a subproblem
with those values fixed and generates either:
- **Optimality cuts** when the subproblem is feasible, or
- **Feasibility cuts** when the subproblem is infeasible and a dual certificate
  is available.

## Construction

```julia
ClassicalOracle(
    data::AbstractData,
    master::Master;
    customize = customize_sub_model!,
    scen_idx::Int = 0,
    param::ClassicalOracleParam = ClassicalOracleParam(),
)
```

### Arguments
- `data`: User-defined problem data passed to the subproblem builder.
- `master`: The master module, whose coupling variables are copied into the
subproblem.
- `customize`: A user-provided function that builds the subproblem model. It is
called as
```julia
customize(model, data, scen_idx; x_copy...)
```
where `x_copy` are the copied master variables.
- `scen_idx`: Scenario index, useful for stochastic or multi-scenario models.
- `param`: Oracle parameters controlling numerical tolerances.

## Fields
- `param`: Oracle parameters (ClassicalOracleParam).
- `model`: The JuMP subproblem model.
- `fixed_x_constraints`: Constraints fixing copied master variables to the
current master solution during cut generation.

## Notes
Master variables are copied into the subproblem with identical axes and names.
Linking constraints of the form master variable copy == value are added and updated at each
oracle call.
Dual values of these constraints are used to construct Benders cuts.
See also: `generate_cuts`, `Master`, `AbstractTypicalOracle`
"""
mutable struct ClassicalOracle <: AbstractTypicalOracle
    
    param::ClassicalOracleParam

    model::Model
    fixed_x_constraints::Vector{ConstraintRef}

    function ClassicalOracle(data::AbstractData, master::Master; 
                            customize = customize_sub_model!,
                            scen_idx::Int=0, 
                            param::ClassicalOracleParam = ClassicalOracleParam())
    
            @debug "Building classical oracle"
            model = Model()

            # Copy the master’s coupling variables into the submodel (with identical axes and symbols)
            x_copy = copy_variables!(model, master.x_tuple)

            # Build the submodel using user-defined customization, passing the copied variables
            customize(model, data, scen_idx; x_copy...)

            # Collect all copied master variables and add linking constraint
            x = var_from_tuple(x_copy)
            @constraint(model, fix_x, x .== 0)

            new(param, model, fix_x)
    end

    ClassicalOracle() = new()
end

function generate_cuts(oracle::ClassicalOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    set_time_limit_sec(oracle.model, time_limit)
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    optimize!(oracle.model)
    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end
    
    status = dual_status(oracle.model)
    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(oracle.model)

        a_x = dual.(oracle.fixed_x_constraints) 
        a_t = [-1.0] 
        a_0 = sub_obj_val - a_x'*x_value 
        if sub_obj_val >= t_value[1] * (1 + oracle.param.rtol) + oracle.param.atol / tol_normalize
            return false, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        else
            return true, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        end

    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(oracle.model)
            dual_sub_obj_val = dual_objective_value(oracle.model)
            @debug "dual_sub_obj_val = $dual_sub_obj_val"
            a_x = dual.(oracle.fixed_x_constraints)
            a_t = [0.0]
            a_0 = dual_sub_obj_val - a_x' * x_value 
            if dual_sub_obj_val >= oracle.param.zero_tol / tol_normalize
                return false, [Hyperplane(a_x, a_t, a_0)], [Inf]
            else
                return true, [Hyperplane(a_x, a_t, a_0)], [Inf]
            end
        end
    else
        throw(UnexpectedModelStatusException("ClassicalOracle: $(status)"))
    end
end




