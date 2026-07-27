"""
    transfer_scaled_linear_rows_and_bounds_with_types!(master, x, dcglp, omega, omega0)

Copy linear rows and scalar bounds from a master model into the DCGLP
formulation used by disjunctive DCGLP oracles.

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
