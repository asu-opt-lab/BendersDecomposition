export customize_master_model!, customize_sub_model!, copy_variables!, var_from_tuple

"""
    customize_master_model!(model::Model, problem::AbstractData)

User-defined hook for constructing the master model.

This function is intended to be implemented by users for their specific
problem type. Given a JuMP `model` and a problem instance data `problem`,
the method should:

1. Add all master-level variables to the model.
2. Add all master-level constraints and any auxiliary structures needed.
3. Return a `NamedTuple` containing the master decision variables  
   (for use by the oracles) and a vector of any auxiliary variables needed by the
   algorithm.

By default, this function throws an `UndefError`, indicating that no
implementation exists for the given subtype of `AbstractData`. Users must
provide a specialized method:

```julia
function customize_master_model!(model::Model, problem::MyDataType)
    # build variables and constraints
    @variable(model, u[1:10])
    @variable(model, v[1:2, 2:3, [:A,:B], 4:5])
    @variable(model, t[1:10]) # auxiliary variables for Benders

    # return variables in a NamedTuple
    return (u = u, v = v), t
end

Arguments
- `model::Model`: A JuMP model that will be modified in place.
- `problem::AbstractData`: A user-defined data object describing the instance required to build the formulation.

Returns
- A `NamedTuple` mapping variable symbolic names to the JuMP variable containers (e.g., `Vector{VariableRef}`, `DenseAxisArray`, or `SparseAxisArray`) created for the master problem.

Notes
- This function must be implemented by the user for each concrete
subtype of `AbstractData`.

"""
function customize_master_model!(model::Model, problem::AbstractData)
    throw(UndefError("update customize_master_model! for $(typeof(problem))"))
end

"""
    customize_sub_model!(model::Model, problem::AbstractData, scen_idx::Int; kwargs...)

User-defined hook for constructing a subproblem model associated with a
specific scenario.

This function must be implemented by users for their concrete subtype of
`AbstractData`. Given a JuMP `model`, a problem instance `problem`, and a
scenario index `scen_idx`, the method should:

1. Add all subproblem variables associated with the given scenario.
2. Add all subproblem constraints associated with the given scenario,
   using master variables passed through `kwargs` (intended for oracle use).

By default, this method throws an `UndefError`, indicating that no implementation
exists for the provided argument types. Users must define their own specialized
method. For example:

```julia
function customize_sub_model!(model::Model, problem::MyDataType, scen_idx::Int; u, v)
    # Create variables
    @variable(model, y >= 0)

    scen_data = problem[scen_idx]

    # Add constraints based on scenario data, referencing master variables u and v
    @constraint(model, sum(u) + y == 1)
    @constraint(model, sum(v) == 1)
    @constraint(model, v[2, 3, :A, 4] + y == 1)
end

Arguments
- `model::Model` A JuMP model that will be modified in place to represent the subproblem.
- `problem::AbstractData` A user-defined data object describing the information required to build the
subproblem formulation.
- `scen_idx::Int` Index of the scenario for which the subproblem is built. This may be -1 if the model does not explicitly use scenarios.
- `kwargs...` Symbolic names of the master variables passed from the master model to the subproblem for use by the oracle.

Notes
- This function must be implemented by the user to construct problem-specific sub-models.
"""

function customize_sub_model!(model::Model, problem::AbstractData, scen_idx::Int; kwargs...) 
    throw(UndefError("update customize_sub_model! for $(typeof(problem))"))
end
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
