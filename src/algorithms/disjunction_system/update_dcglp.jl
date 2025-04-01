function get_subproblem_value(sub::AbstractSubProblem, x_value::Vector{Float64}, cut_strategy::CutStrategy)
    
    # # Check if master.x_value is close enough to integer values
    if !all(x -> isapprox(x, round(x), atol=1e-4), x_value)
        return Inf
    end
    
    if dual_status(sub.model) == FEASIBLE_POINT
        return objective_value(sub.model)
    elseif dual_status(sub.model) == INFEASIBLE_POINT
        return Inf
    else
        error("Subproblem is not feasible or optimal")
    end
end

function get_subproblem_value(sub::KnapsackUFLPSubProblem, x_value::Vector{Float64}, cut_strategy::FatKnapsackCut)
    _, sub_obj_val = generate_cut_coefficients(sub, x_value, cut_strategy)
    return sub_obj_val  
end

function select_disjunctive_inequality(x_value)

    gap_x = @. abs(x_value - 0.5)
    index = argmin(gap_x)
    a = zeros(Int, length(x_value))
    a[index] = 1
    @info "Selected disjunction index: $index, gap_x: $(x_value[index])"
    return a, 0
end

function update_dcglp!(dcglp::DCGLP, disjunctive_inequality::Tuple{Vector{Int}, Int}, disjunction_system::DisjunctiveCut)
    replace_disjunctive_inequality!(dcglp, disjunctive_inequality, disjunction_system.norm_type)
    update_added_benders_constraints!(dcglp, disjunction_system)
    # add_disjunctive_cut!(dcglp)
end


function replace_disjunctive_inequality!(dcglp::DCGLP, disjunctive_inequality, ::LNorm)
    # Delete previous constraints if they exist
    if !isempty(dcglp.disjunctive_inequalities_constraints)
        for constraint in dcglp.disjunctive_inequalities_constraints
            delete(dcglp.model, constraint)
        end
        empty!(dcglp.disjunctive_inequalities_constraints)
    end

    # Add new constraints
    coef, constant = disjunctive_inequality
    push!(dcglp.disjunctive_inequalities_constraints, 
        @constraint(dcglp.model, 0 >= dcglp.model[:k₀]*(constant+1) - coef'dcglp.model[:kₓ]))
    push!(dcglp.disjunctive_inequalities_constraints,
        @constraint(dcglp.model, 0 >= -dcglp.model[:v₀]*constant + coef'dcglp.model[:vₓ]))
end

function update_added_benders_constraints!(dcglp::DCGLP, disjunctive_system::DisjunctiveCut)
    if !disjunctive_system.reuse_dcglp
        for constraint in dcglp.dcglp_constraints
            delete(dcglp.model, constraint)
        end
    end
    empty!(dcglp.dcglp_constraints)
    empty!(dcglp.master_cuts)
end


function add_disjunctive_cut!(dcglp::DCGLP)
    if !isempty(dcglp.γ_values)
        γ₀, γₓ, γₜ = dcglp.γ_values[end]
        @constraint(dcglp.model, γ₀*dcglp.model[:k₀] + dot(γₓ, dcglp.model[:kₓ]) + dot(γₜ, dcglp.model[:kₜ]) <= 0)
        @constraint(dcglp.model, γ₀*dcglp.model[:v₀] + dot(γₓ, dcglp.model[:vₓ]) + dot(γₜ, dcglp.model[:vₜ]) <= 0)    
    end
end



# for testing

# function update_dcglp!(algo::AbstractBendersAlgorithm{D,L,DisjunctiveCut}, disjunctive_inequality::Tuple{Vector{Int}, Int}) where {D<:AbstractData, L<:SequentialLoop}
#     algo.dcglp = create_dcglp(algo.data, algo.cut_strategy)
#     set_optimizer(algo.dcglp.model, CPLEX.Optimizer)
#     set_silent(algo.dcglp.model)
#     coef, constant = disjunctive_inequality
#     @constraint(algo.dcglp.model, 0 >= algo.dcglp.model[:k₀]*(constant+1) - coef'algo.dcglp.model[:kₓ])
#     @constraint(algo.dcglp.model, 0 >= -algo.dcglp.model[:v₀]*constant + coef'algo.dcglp.model[:vₓ])
#     add_disjunctive_cuts!(algo.dcglp)
# end

# function add_disjunctive_cuts!(dcglp::DCGLP)
#     for γ in dcglp.γ_values

#         γ₀, γₓ, γₜ = γ
#         @constraint(dcglp.model, γ₀*dcglp.model[:k₀] + sum(γₓ.*dcglp.model[:kₓ]) + γₜ*dcglp.model[:kₜ] <= 0)
#         @constraint(dcglp.model, γ₀*dcglp.model[:v₀] + sum(γₓ.*dcglp.model[:vₓ]) + γₜ*dcglp.model[:vₜ] <= 0)
#     end
# end