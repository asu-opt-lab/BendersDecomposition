export transfer_scaled_linear_rows_and_bounds_with_types!


function transfer_scaled_linear_rows_and_bounds_with_types!(
    master::Model,
    x::Vector{VariableRef},
    dcglp::Model,
    omega::Vector{VariableRef},
    omega0::VariableRef
)
    # detect unhandled constraint types and issue warnings
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
        @warn "A master constraint of type ($F, $S) was not automatically incorporated into dcglp. \
        If this constraint is linear, please add it manually."
    end

    # Build mapping from MOI variable index -> position in x/omega
    idx_to_pos = Dict{Int,Int}()
    for (pos, v) in enumerate(x)
        vi = JuMP.index(v)            # MOI.VariableIndex
        idx_to_pos[vi.value] = pos
    end
    
    # basic sanity
    if length(x) != length(omega)
        error("x and omega must have the same length/structure.")
    end

    backend = JuMP.backend(master)

    pair_types = [
        (MOI.VariableIndex, MOI.GreaterThan{Float64}),
        (MOI.VariableIndex, MOI.LessThan{Float64}),
        (MOI.VariableIndex, MOI.EqualTo{Float64}),
        (MOI.VariableIndex, MOI.Interval{Float64}),
        (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}),
        (MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}),
        (MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}),
        (MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64})
    ]

    for (F,S) in pair_types
        cis = MOI.get(backend, MOI.ListOfConstraintIndices{F,S}())
        for ci in cis
            f = MOI.get(backend, MOI.ConstraintFunction(), ci) # ScalarAffineFunction
            s = MOI.get(backend, MOI.ConstraintSet(), ci)      # set type
            # build LHS as AffExpr (sum over terms mapped to omega)
            if f isa MOI.ScalarAffineFunction
                if isempty(f.terms)
                    lhs = 0.0
                else
                    # accumulate as AffExpr: start from zero and add terms
                    lhs = zero(AffExpr)
                    for term in f.terms
                        vi_val = term.variable.value
                        pos = get(idx_to_pos, vi_val, nothing)
                        if pos === nothing
                            break # skip if the constraint contains a variable not present in `x`
                        end
                        coeff = term.coefficient
                        lhs += coeff * omega[pos]
                    end
                end
            elseif f isa MOI.VariableIndex
                pos = get(idx_to_pos, f.value, nothing)
                if pos === nothing
                    break # skip if the constraint contains a variable not present in `x`
                end
                lhs = omega[pos]
            end

            added_constraints = []
            # determine bounds present for this set type
            if S === MOI.GreaterThan{Float64}
                lower = s.lower
                # add lhs >= lower * omega0
                push!(added_constraints, @constraint(dcglp, lhs >= lower * omega0))
            elseif S === MOI.LessThan{Float64}
                upper = s.upper
                push!(added_constraints, @constraint(dcglp, lhs <= upper * omega0))
            elseif S === MOI.EqualTo{Float64}
                val = s.value
                push!(added_constraints, @constraint(dcglp, lhs == val * omega0))
            elseif S === MOI.Interval{Float64}
                lower = s.lower
                upper = s.upper
                if lower == upper
                    push!(added_constraints, @constraint(dcglp, lhs == lower * omega0))
                else
                    push!(added_constraints, @constraint(dcglp, lhs >= lower * omega0))
                    push!(added_constraints, @constraint(dcglp, lhs <= upper * omega0))
                end
            end
            @show added_constraints
        end
    end
end