"""
    _validate_lp_compatibility(model::Model)

Validate that the model is compatible with typical oracles.
Typical oracles require a Linear Programming (LP) subproblem:
- no discontinuous (e.g., integer/binary/semi-*) variables,
- affine objective and constraints.

Throws `UnsupportedModelException` when unsupported model components are found.

This validation is shared by [`ClassicalOracle`](@ref), [`UnifiedOracle`](@ref), and [`ParetoOracle`](@ref).
"""
function _validate_lp_compatibility(model::Model)
    # 1. Check for discontinuous variables
    for (_, S) in list_of_constraint_types(model)
        if S <: Union{MOI.Integer, MOI.ZeroOne, MOI.Semicontinuous, MOI.Semiinteger}
            throw(UnsupportedModelException(
                "_validate_lp_compatibility: Unsupported variable type: $S. " *
                "Typical oracles require a Linear Programming (LP) subproblem. " *
                "Discontinuous variables are not allowed."
            ))
        end
    end

    # 2. Check for non-affine objective function
    obj_type = objective_function_type(model)
    if !(obj_type <: Union{VariableRef, AffExpr, Real})
        throw(UnsupportedModelException(
            "_validate_lp_compatibility: Unsupported objective function type: $obj_type. " *
            "Typical oracles require a Linear Programming (LP) subproblem. "
        ))
    end

    # 3. Check for non-affine constraints
    for con in all_constraints(model, include_variable_in_set_constraints=true)
        con_obj = constraint_object(con)
        set = con_obj.set
        func = con_obj.func

        if func isa VariableRef || func isa AffExpr || (func isa AbstractVector && all(x -> x isa AffExpr || x isa VariableRef, func))
            if !(set isa MOI.GreaterThan || set isa MOI.LessThan || set isa MOI.EqualTo || set isa MOI.Interval || set isa MOI.Zeros || set isa MOI.Nonpositives || set isa MOI.Nonnegatives)
                throw(UnsupportedModelException(
                    "_validate_lp_compatibility: Unsupported constraint set type: $(typeof(set)). " *
                    "Typical oracles only support affine constraints."
                ))
            end
        else
            throw(UnsupportedModelException(
                "_validate_lp_compatibility: Unsupported constraint function type: $(typeof(func)). " *
                "Typical oracles only support affine constraints."
            ))
        end
    end
end

"""
    copy_variables!(model::Model, x::NamedTuple) -> NamedTuple

Create JuMP variables inside `model` that mirror the structure of the NamedTuple `x`.

For each `(key, value)` pair in `x` the function:

- If `value` is a `JuMP.VariableRef`, creates a new scalar variable and stores it at `model[key]`.
- If `value` is a `JuMP.Containers.Array` or `JuMP.Containers.DenseAxisArray`, creates a JuMP array variable with the same axes and stores it at `model[key]`.
- If `value` is a `JuMP.Containers.SparseAxisArray`, creates a sparse variable container with the same index keys and stores it at `model[key]`.
- Otherwise logs an error (`@error`) and does not add a variable for that key.

# Arguments
- `model::Model` : a JuMP model that will be mutated (new variables are added).
- `x::NamedTuple` : named tuple whose values are JuMP variable containers or `VariableRef`s from another JuMP model.

# Returns
A `NamedTuple` with the same keys as `x` whose values are the newly-created JuMP variable(s) for `model`.
The function registers each created variable container on `model` under the original `key`.

# Example
```julia
using JuMP, JuMP.Containers

m = Model()
# original variables (in a different model)
@variable(m, a >= 0)                    # scalar VariableRef
@variable(m, b[2:3])                             # DenseAxisArray
@variable(m, c[i=1:2, j=i:2, k=j:4])             # SparseAxisArray
x = (a = a, b = b, c = c)

model2 = Model()

# Create variables in `model2` mirroring structure of `x`
xvars = copy_variables!(model2, x)

# Now `xvars.a`, `xvars.b`, `xvars.c` are JuMP variables created inside `model2`.
"""
function copy_variables!(model::Model, x::NamedTuple)
    for (key, value) in pairs(x)
        if value isa JuMP.VariableRef
            model[key] = @variable(model, base_name = "$(key)")
        elseif value isa JuMP.Containers.Array
            # axes returns a tuple of ranges, e.g. (1:3, 2:4) for x[1:3, 2:4]
            arr = Array{VariableRef}(undef, size(value))

            for idx in keys(value)
                i = idx isa Number ? idx : idx.I
                arr[idx] = @variable(model, base_name = "$(key)[$(i)]")
            end
            model[key] = arr
        elseif value isa JuMP.Containers.DenseAxisArray
            # ranges = [axes(value)...]
            # expr = :( @variable($(model), $(key)[$(ranges...)], base_name = $(string(key))) )
            # eval(expr)

            # axes returns a tuple of ranges, e.g. (1:3, 2:4) for x[1:3, 2:4]
            ax = axes(value)
            vartemp = Array{VariableRef}(undef, size(value))
            darr = Containers.DenseAxisArray(vartemp, ax...)

            for idx in keys(value)
                i = idx isa Number ? idx : idx.I
                darr[idx] = @variable(model, base_name = "$(key)[$(i)]")
            end
            model[key] = darr
        elseif value isa JuMP.Containers.SparseAxisArray
            vartemp = @variable(model, [eachindex(value)], base_name = "$(key)")
            varmap = Dict(index => vartemp[index] for index in eachindex(value))
            var = Containers.SparseAxisArray(varmap)
            model[key] = var
        else
            @error "We are not currently accepting a custom JuMP variable container: $(typeof(value))"
        end
    end

    x_copy = (; (key => model[key] for (key, _) in pairs(x))...)

    return x_copy
end

"""
    var_from_tuple(x_tuple::NamedTuple) -> Vector{VariableRef}

Extract all JuMP variables contained in a `NamedTuple`—including scalars, Arrays, `DenseAxisArray`s, and `SparseAxisArray`s—and return them as a flat `Vector{VariableRef}`.
"""
function var_from_tuple(x_tuple::NamedTuple)
    x = Vector{VariableRef}(undef,0)
    for value in values(x_tuple)
        if value isa Array
            # for Array Container
            append!(x, collect(value))
        else
            if value.data isa Array
                # for DenseAxisArray Container
                append!(x, collect(value.data))
            else
                # for SparseAxisArray Container
                append!(x, collect(values(value.data)))
            # else
            #     @error "unexpected Variable Container values"
            end
        end
    end
    return x
end
