"""
    update_master_model!(model::Model, data::AbstractData) -> (x, t)

Build one master problem formulation in `model`.

This function tells BendersX how to build the master problem for your data type.
Define one method for your custom data type:

```julia
function update_master_model!(model::Model, data::MyDataType)
    # Add variables, constraints, and objective.
    return x, t
end
```

The function receives an empty JuMP `model`. It should:

1. Create the master variables.
2. Add the master constraints and objective.
3. Return `(x, t)`.

`x` must be a `NamedTuple` containing the master variables that subproblems may
need to read. The field names in `x` are the names subproblems use.

`t` is the auxiliary variable, or collection of variables, that BendersX uses
when adding Benders cuts. The variable does not have to be named `t` in your
model; only the returned value matters.

Optimizer selection is handled elsewhere. Do not attach an optimizer in this
function.

# Example

```julia
struct MyDataType <: AbstractData
    n::Int
end

function update_master_model!(model::Model, data::MyDataType)
    @variable(model, u[1:data.n])
    @variable(model, theta)

    @constraint(model, sum(u) >= 1)
    @objective(model, Min, theta)

    return (u = u,), theta
end
```

You may also pass a custom builder directly:

```julia
Master(data; model = build_master_model!)
```

In that case, `build_master_model!(model, data)` must follow the same rule:
modify `model` and return `(x, t)`.
"""
function update_master_model!(model::Model, data::AbstractData)
    throw(UndefError(
        "BendersX does not know how to build a master model for $(typeof(data)). " *
        "Define `update_master_model!(model::Model, data::$(typeof(data)))` so " *
        "`Master(data)` can use it as the default builder, or pass a builder " *
        "explicitly with `Master(data; model = build_master_model!)`. The builder " *
        "can have any name, but it must accept `(model, data)` and return `(x, t)`."
    ))
end

"""
    update_sub_model!(model::Model, data::AbstractData, scen_idx::Int; kwargs...)

Build one subproblem formulation in `model`.

This function tells BendersX how to build a subproblem for your data type.
Define one method for your custom data type:

```julia
function update_sub_model!(model::Model, data::MyDataType, scen_idx::Int; x)
    # Add variables, constraints, and objective.
    return nothing
end
```

The function receives an empty JuMP `model`. It should:

1. Create the subproblem variables.
2. Add the subproblem constraints and objective.
3. Use the master variables passed as keyword arguments when the subproblem
   depends on the master solution.

`scen_idx` identifies which scenario is being built. If your model is not
scenario-based, you can ignore this argument.

The keyword arguments are the master variables returned by
`update_master_model!`. For example, if the master returns `(x = x,)`, then the
subproblem can accept `; x` and use `x` directly in JuMP expressions.

Most subproblem builders return `nothing`. If you need generalized bound
constraints, return `(gbc_lhs, gbc_rhs, gbc_sense)`, where each relation has the
form `subproblem_variable <=/>=/== affine_expression_of_master_variables`.

Optimizer selection is handled elsewhere. Do not attach an optimizer in this
function.

# Example

```julia
struct MyDataType <: AbstractData
    demand::Vector{Float64}
    cost::Vector{Float64}
end

function update_sub_model!(model::Model, data::MyDataType, scen_idx::Int; x)
    @variable(model, y[eachindex(data.demand)] >= 0)

    @constraint(model, demand[j in eachindex(data.demand)], y[j] == data.demand[j])
    @constraint(model, linking[j in eachindex(data.demand)], y[j] <= x[j])
    @objective(model, Min, sum(data.cost[j] * y[j] for j in eachindex(data.demand)))

    return nothing
end
```

You may also pass a custom builder directly:

```julia
ClassicalOracle(data, master; model = build_sub_model!)
```

In that case, `build_sub_model!(model, data, scen_idx; kwargs...)` must follow
the same rule: modify `model` and return `nothing`, or return generalized bound
constraint data when needed.
"""
function update_sub_model!(model::Model, data::AbstractData, scen_idx::Int; kwargs...)
    throw(UndefError(
        "BendersX does not know how to build a subproblem model for $(typeof(data)). " *
        "Define `update_sub_model!(model::Model, data::$(typeof(data)), scen_idx::Int; kwargs...)` " *
        "to construct a model-based oracle without `model = ...`, or pass your own " *
        "subproblem-model function with `model = build_sub_model!` when constructing the oracle."
    ))
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

"""
    transfer_scaled_linear_rows_and_bounds_with_types!(master, x, dcglp, omega, omega0)

Copy linear rows and scalar bounds from a master model into the DCGLP
formulation used by [`SplitOracle`](@ref).

Only linear constraints whose variables all belong to the supplied vector `x`
are transferred. Variable integrality restrictions are ignored, because the
target DCGLP is a continuous lifting of the master model.
"""
function transfer_scaled_linear_rows_and_bounds_with_types!(
    master::Model,
    x::Vector{VariableRef},
    dcglp::Model,
    omega::Vector{VariableRef},
    omega0::VariableRef,
)
    pairs_present = JuMP.list_of_constraint_types(master)
    for (F, S) in pairs_present
        if F in [AffExpr; VariableRef]
            if S in [MOI.GreaterThan{Float64}; MOI.LessThan{Float64}; MOI.EqualTo{Float64}; MOI.Interval{Float64}]
                continue
            end
        end
        if S in [MOI.Integer; MOI.ZeroOne]
            continue
        end
        @warn "A master constraint of type ($F, $S) was not automatically incorporated into dcglp. If this constraint is linear, please add it manually."
    end

    idx_to_pos = Dict{Int,Int}()
    for (pos, v) in enumerate(x)
        vi = JuMP.index(v)
        idx_to_pos[vi.value] = pos
    end

    length(x) == length(omega) || error("x and omega must have the same length/structure.")

    backend = JuMP.backend(master)
    pair_types = [
        (MOI.VariableIndex, MOI.GreaterThan{Float64}),
        (MOI.VariableIndex, MOI.LessThan{Float64}),
        (MOI.VariableIndex, MOI.EqualTo{Float64}),
        (MOI.VariableIndex, MOI.Interval{Float64}),
        (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}),
        (MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}),
        (MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}),
        (MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}),
    ]

    for (F, S) in pair_types
        for ci in MOI.get(backend, MOI.ListOfConstraintIndices{F,S}())
            func = MOI.get(backend, MOI.ConstraintFunction(), ci)
            set = MOI.get(backend, MOI.ConstraintSet(), ci)

            terms = Tuple{Float64,Int}[]
            constant = 0.0

            if F == MOI.VariableIndex
                vpos = get(idx_to_pos, func.value, 0)
                vpos == 0 && continue
                push!(terms, (1.0, vpos))
            else
                constant = func.constant
                ok = true
                for term in func.terms
                    vpos = get(idx_to_pos, term.variable.value, 0)
                    if vpos == 0
                        ok = false
                        break
                    end
                    push!(terms, (term.coefficient, vpos))
                end
                ok || continue
            end

            expr = sum(a * omega[j] for (a, j) in terms)

            # MOI stores scalar affine constraints as `constant + expr in set`.
            # The lifted point (omega / omega0) must satisfy the same relation.
            if S == MOI.GreaterThan{Float64}
                @constraint(dcglp, expr + constant * omega0 >= set.lower * omega0)
            elseif S == MOI.LessThan{Float64}
                @constraint(dcglp, expr + constant * omega0 <= set.upper * omega0)
            elseif S == MOI.EqualTo{Float64}
                @constraint(dcglp, expr + constant * omega0 == set.value * omega0)
            elseif S == MOI.Interval{Float64}
                @constraint(dcglp, expr + constant * omega0 <= set.upper * omega0)
                @constraint(dcglp, expr + constant * omega0 >= set.lower * omega0)
            end
        end
    end
end
