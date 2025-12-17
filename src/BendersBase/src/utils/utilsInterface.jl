export customize_master_model!, customize_sub_model!, copy_variables!, var_from_tuple

"""
    customize_master_model!(model::Model, data::AbstractData) -> NamedTuple, Vector{VariableRef}

User-defined hook for constructing the master model.

This function is intended to be implemented by users for their specific
problem. Given a JuMP `model` and an instance `data`,
the method should:

1. Add all master-level variables to the model.
2. Return a `NamedTuple` containing the master decision variables x 
and a vector of any auxiliary variables t (for use by the oracles).

By default, this function throws an `UndefError`, indicating that no
implementation exists for the given subtype of `AbstractData`. Users must
provide a specialized method:

```julia
struct MyDataType <: AbstractData
    ...
end

function customize_master_model!(model::Model, data::MyDataType)
    # build variables and constraints
    @variable(model, u[1:10])
    @variable(model, v[1:2, 2:3, [:A,:B], 4:5])
    @variable(model, t[1:10]) # auxiliary variables for Benders

    # return x-variables in a NamedTuple and t-variables in a vector
    return (u = u, v = v), t
end

Arguments
- `model::Model`: A JuMP model that will be modified in place.
- `data::MyDataType`: A user-defined data object describing the instance required to build the formulation.

Returns
- A `NamedTuple` mapping variable symbolic names to the JuMP variable containers (e.g., `Vector{VariableRef}`, `DenseAxisArray`, or `SparseAxisArray`) created for the master problem.

Notes
- This function must be implemented by the user for each concrete
subtype of `AbstractData`.

"""
function customize_master_model!(model::Model, data::AbstractData)
    throw(UndefError("update customize_master_model! for $(typeof(data))"))
end

"""
    customize_sub_model!(model::Model, data::AbstractData, scen_idx::Int; kwargs...)

User-defined hook for constructing a subproblem model associated with a
specific scenario.

This function must be implemented by users for their concrete subtype of
`AbstractData`. Given a JuMP `model`, a problem instance `data`, and a
scenario index `scen_idx`, the method should:

1. Add all subproblem variables associated with the given scenario.
2. Add all subproblem constraints associated with the given scenario,
   using master variables passed through `kwargs`.

By default, this method throws an `UndefError`, indicating that no implementation
exists for the provided argument types. Users must define their own specialized
method if they want to use any model-based oracle. For example:

```julia
function customize_sub_model!(model::Model, data::MyDataType, scen_idx::Int; u, v)
    # Create variables
    @variable(model, y >= 0)

    scen_data = data[scen_idx]

    # Add constraints based on scenario data, referencing master variables u and v
    @constraint(model, sum(u) + y == 1)
    @constraint(model, sum(v) == 1)
    @constraint(model, v[2, 3, :A, 4] + y == 1)
end

Arguments
- `model::Model` A JuMP model that will be modified in place to represent the subproblem.
- `data::AbstractData` A user-defined data object describing the information required to build the
subproblem formulation.
- `scen_idx::Int` Index of the scenario for which the subproblem is built. This may be -1 if the model does not explicitly use scenarios.
- `kwargs...` Symbolic names of the master variables passed from the master model to the subproblem for use by the oracle.

Notes
- This function must be implemented by the user to construct model-based oracles.
"""
function customize_sub_model!(model::Model, data::AbstractData, scen_idx::Int; kwargs...)
    throw(UndefError("update customize_sub_model! for $(typeof(data))"))
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
