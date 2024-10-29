export DisjunctionSystem

mutable struct DisjunctionSystem <: CutGenerationStrategy
    norm_type::AbstractNormType
    _cut_strategy::CutGenerationStrategy
end

function generate_cuts(algo::AbstractBendersAlgorithm, disjunction_system::DisjunctionSystem)
    disjunctive_inequality = select_disjunctive_inequality(algo.master.x_value)
    dcglp_problem = create_dcglp(algo.data, disjunctive_inequality, disjunction_system._cut_strategy, disjunction_system.norm_type)
    assign_solver!(dcglp_problem.model, :CPLEX)
    solve!(dcglp_problem, algo)

    γ₀, γₓ, γₜ = dual(dcglp_problem.γconstraints[:γ₀]), dual.(dcglp_problem.γconstraints[:γₓ]), dual.(dcglp_problem.γconstraints[:γₜ])
    cut = @expression(algo.master.model, γ₀ + γₓ'algo.master.model[:x] + γₜ'algo.master.model[:t])

    _, obj_value = generate_cuts(algo, disjunction_system._cut_strategy)
    return cut, obj_value
end

function solve!(dcglp_problem::AbstractDCGLP, algo::AbstractBendersAlgorithm)
    LB, UB = -Inf, Inf
    LBs = []
    iter = 0
    x_value = algo.master.x_value
    t_value = value.(algo.master.var[:t])
    set_normalized_rhs.(dcglp_problem.model[:conx], x_value)
    set_normalized_rhs.(dcglp_problem.model[:cont], t_value)
    start_time = time()
    while true

        iter += 1
        optimize!(dcglp_problem.model)
        k̂₀ = value(dcglp_problem.model[:k₀])
        k̂ₓ = value.(dcglp_problem.model[:kₓ])
        k̂ₜ = value.(dcglp_problem.model[:kₜ])
        v̂₀ = value(dcglp_problem.model[:v₀])
        v̂ₓ = value.(dcglp_problem.model[:vₓ])
        v̂ₜ = value.(dcglp_problem.model[:vₜ])
        τ̂ = value(dcglp_problem.model[:τ])
        _sx = value.(dcglp_problem.model[:sx])

        obj_value_k, obj_value_v = Inf, Inf
        if k̂₀ != 0 
            dual_values_k, obj_value_k = generate_cut_coefficients(algo.sub, k̂ₓ./k̂₀, algo.cut_strategy._cut_strategy)
        end

        if v̂₀ != 0 
            dual_values_v, obj_value_v = generate_cut_coefficients(algo.sub, v̂ₓ./v̂₀, algo.cut_strategy._cut_strategy)
        end


        ##################### LB and UB #####################
        LB = τ̂
        UB = update_UB!(UB,_sx,obj_value_k,obj_value_v,t_value, algo.cut_strategy.norm_type)
        push!(LBs, LB)

        @info "Iteration $iter: LB = $LB, UB = $UB, _UB1 = $obj_value_k, _UB2 = $obj_value_v"


        ##################### check termination #####################
        if (UB - LB)/abs(UB) <= 1e-2 || (UB - LB) <= 0.01
            @info "Optimal solution found"
            break
        end

        if iter >= 5 && all(LBs[end] - LB <= 1e-05 for LB in LBs[end-4:end-1]) || 200 <= time() - start_time
            @info "Time limit reached"
            break
        end

        ##################### add cuts into DCGLP and master problem #####################
        if k̂₀ != 0
            if obj_value_k >= 1e-3
                cuts_k, cuts_v, cuts_master = build_cuts(dcglp_problem, algo.master, algo.sub, dual_values_k, algo.cut_strategy._cut_strategy)
                @constraint(dcglp_problem.model, 0 .>= cuts_k)
                @constraint(dcglp_problem.model, 0 .>= cuts_v)
                @constraint(algo.master.model, 0 .>= cuts_master)
            end
        end

        if v̂₀ != 0
            if obj_value_v >= 1e-3
                cuts_k, cuts_v, cuts_master = build_cuts(dcglp_problem, algo.master, algo.sub, dual_values_v, algo.cut_strategy._cut_strategy)
                @constraint(dcglp_problem.model, 0 .>= cuts_k)
                @constraint(dcglp_problem.model, 0 .>= cuts_v)
                @constraint(algo.master.model, 0 .>= cuts_master)
            end
        end
        
    end

    ##################### update #####################

end

function update_UB!(UB,_sx,obj_value_k,obj_value_v,t_value,::L1Norm) return min(UB,norm([ _sx; obj_value_k .+ obj_value_v .- t_value], Inf)) end
function update_UB!(UB,_sx,obj_value_k,obj_value_v,t_value,::L2Norm) return min(UB,norm([ _sx; obj_value_k .+ obj_value_v .- t_value], 2)) end
function update_UB!(UB,_sx,obj_value_k,obj_value_v,t_value,::LInfNorm) return min(UB,norm([ _sx; obj_value_k .+ obj_value_v .- t_value], 1)) end


function select_disjunctive_inequality(x_value)
    # Calculate the gap between each x value and 0.5
    gap_x = abs.(x_value .- 0.5)
    
    # Find the index of the x value closest to 0.5
    index = argmin(gap_x)
    
    # Create a vector 'a' with zeros, except for a 1 at the found index
    a = zeros(Int, length(x_value))
    a[index] = 1
    
    # Log the selected index for debugging
    @debug "Selected disjunction index: $index"
    
    # Return the disjunction vector 'a' and the right-hand side 0
    return a, 0
end

