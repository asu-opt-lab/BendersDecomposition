"""
    _insert_suffix(name::String, suffix::String) -> String

Insert suffix before the index brackets in a constraint name.
Examples:
- `demand[1]` + `_lb` -> `demand_lb[1]`
- `flow[i,j]` + `_ub` -> `flow_ub[i,j]`
- `simple` + `_lb` -> `simple_lb`
"""
function _insert_suffix(name::String, suffix::String)
    bracket_pos = findfirst('[', name)
    if bracket_pos === nothing
        return name * suffix
    else
        return name[1:bracket_pos-1] * suffix * name[bracket_pos:end]
    end
end

"""
    _scalarize_constraints!(model::Model)

Convert all non-scalar (e.g., vectorized) constraints in `model` into equivalent
scalar constraints so the resulting model contains only MOI constraint
types `GreaterThan`, `LessThan`, and `EqualTo`.

Splits:
- `MOI.Interval` → `GreaterThan(lower)` + `LessThan(upper)`
- `MOI.Zeros` → individual `EqualTo(0)` per row
- `MOI.Nonnegatives` → individual `GreaterThan(0)` per row
- `MOI.Nonpositives` → individual `LessThan(0)` per row

# Example
```julia
    # before: a single vector constraint `A * x >= b`
    _scalarize_constraints!(model)
    # after: N scalar constraints `A[i,:] * x >= b[i]` for i = 1:N
```
"""
function _scalarize_constraints!(model::Model)
    # Collect constraints to split (avoid modifying during iteration)
    cons_to_split = ConstraintRef[]
    for con in all_constraints(model, include_variable_in_set_constraints=true)
        set = constraint_object(con).set
        if set isa MOI.Interval || set isa MOI.Zeros || set isa MOI.Nonnegatives || set isa MOI.Nonpositives
            push!(cons_to_split, con)
        end
    end

    for con in cons_to_split
        con_obj = constraint_object(con)
        set = con_obj.set
        func = con_obj.func
        original_name = name(con)

        if set isa MOI.Interval
            # lb <= expr <= ub  →  expr >= lb  AND  expr <= ub
            lb_con = @constraint(model, func >= set.lower)
            ub_con = @constraint(model, func <= set.upper)
            if !isempty(original_name)
                set_name(lb_con, _insert_suffix(original_name, "_lb"))
                set_name(ub_con, _insert_suffix(original_name, "_ub"))
            end

        elseif set isa MOI.Zeros
            # [f1; f2; ...] ∈ Zeros  →  f1 == 0, f2 == 0, ...
            for (i, fi) in enumerate(func)
                new_con = @constraint(model, fi == 0)
                if !isempty(original_name)
                    set_name(new_con, _insert_suffix(original_name, "_$i"))
                end
            end

        elseif set isa MOI.Nonnegatives
            # [f1; f2; ...] ∈ Nonneg  →  f1 >= 0, f2 >= 0, ...
            for (i, fi) in enumerate(func)
                new_con = @constraint(model, fi >= 0)
                if !isempty(original_name)
                    set_name(new_con, _insert_suffix(original_name, "_$i"))
                end
            end

        elseif set isa MOI.Nonpositives
            # [f1; f2; ...] ∈ Nonpos  →  f1 <= 0, f2 <= 0, ...
            for (i, fi) in enumerate(func)
                new_con = @constraint(model, fi <= 0)
                if !isempty(original_name)
                    set_name(new_con, _insert_suffix(original_name, "_$i"))
                end
            end
        end

        delete(model, con)
    end
end

"""
Append constraints to a JuMP model using a symbolic name.
If constraints with the name `constr_symbol` exist, append to them. Otherwise, create a new constraint group named `constr_symbol`, enforcing `0 .>= exprs`.
"""
function add_constraints(model::Model, constr_symbol::Symbol, exprs::Vector{AffExpr})
    # add constraints in the form of 0 .>= expr
    if haskey(model, constr_symbol)
        append!(model[constr_symbol], @constraint(model, 0 .>= exprs))
    else
        model[constr_symbol] = @constraint(model, 0 .>= exprs)
    end
end

function delete_registered_constraints!(model::Model, sym::Symbol)
    haskey(model, sym) || return nothing
    registered = model[sym]
    if registered isa AbstractArray
        delete.(Ref(model), registered)
    else
        delete(model, registered)
    end
    unregister(model, sym)
end
