export copy_variables!, var_from_tuple

"""
    copy_variables!(model::Model, x::NamedTuple)

Create JuMP variables inside `model` that mirror the structure of the NamedTuple `x`.

For each `(key, value)` pair in `x` the function:

- If `value` is a `JuMP.VariableRef`, creates a new scalar variable and stores it at `model[key]`.
- If `value` is a `JuMP.Containers.Array` or `JuMP.Containers.DenseAxisArray`, creates a JuMP array variable with the same axes and stores it at `model[key]`.
- If `value` is a `JuMP.Containers.SparseAxisArray`, creates a sparse variable container with the same indices and stores it at `model[key]`.
- Otherwise logs an error (`@error`) and does not add a variable for that key.

# Arguments
- `model::Model` : a JuMP model that will be mutated (new variables are added).
- `x::NamedTuple` : named tuple whose values are JuMP variable containers or `VariableRef`s of another JuMP Model.

# Returns
A `NamedTuple` with the same keys as `x` whose values are the newly created JuMP variable(s) for the given model. The function register each created variable container for the `model` with the original key.

# Side effects / Notes
- For dense arrays the implementation builds and `eval`s an expression that calls `@variable` with the original axes. Because `eval` is used, take care when calling this from inside modules or non-global scopes — you may need to adjust the implementation for those contexts.
- For `SparseAxisArray` values the function creates a temporary vector of variables and then maps those to a `Containers.SparseAxisArray`.
- Unsupported JuMP container types will produce an `@error` message; no variable will be created for that key.

# Example
```julia
using JuMP, JuMP.Containers

m = Model()
x = (a = @variable(Model(), base_name="a"),   # scalar VariableRef (example)
     b = Containers.DenseAxisArray{Float64}((1:3,)),  # placeholder array-like shape
     c = Containers.SparseAxisArray(Dict(1=>0.0, 10=>0.0)))  # sparse indices example

# Create variables in `m` mirroring structure of `x`
xvars = copy_variables!(m, x)

# Now `xvars.a`, `xvars.b`, `xvars.c` are JuMP variables (and `m` contains them)
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
