export update_λ, update_y, retrieve_dual_values, gen_figure

"""
utility functions for Dual Decomposition
"""
function update_λ(λ_k, y_k, u, d, x, α)
    @threads for i in eachindex(λ_k)
        y_i = @view y_k[i,:]; subgrad = dot(d, y_i) - u[i]*x[i]
        val = λ_k[i] + α*subgrad
        λ_k[i] = max(0.0, val)
    end
    return λ_k
end

function update_y(y_k, λ_k, c, d, x)
    J = length(d) 
    critical_facility_indices = Vector{Int}(undef, J)
    sorted_indices, c_sorted, x_sorted = sort_costs_indices(J, c + λ_k .* d', x) # sorted_indice(x index): sort x with cheapest order, x_sorted (x value sorted)

    @threads for j in 1:J
        k = find_critical_item(c_sorted[j], x_sorted[j]) # k: index of critical item, sorted_indices[j][k] is a critical item
        critical_facility_indices[j] = k

        jth_sortin = sorted_indices[j]; idxs_prev_k = @view jth_sortin[1:k-1]; idx_after_k = @view jth_sortin[k+1:end]; y_col = @view y_k[:,j]
        
        if k!= 1
            y_col[idxs_prev_k] .= @view x_sorted[j][1:k-1]
            y_col[jth_sortin[k]] = 1 - sum(@view x_sorted[j][1:k-1])
        else
            y_col[jth_sortin[k]] = x_sorted[j][k]
        end
        y_col[idx_after_k] .= 0 
    end
    return y_k, sorted_indices, c_sorted, x_sorted, critical_facility_indices
end

function sort_costs_indices(J, sub_obj, x)
    sorted_indices = Vector{Vector{Int64}}(undef, J)
    c_sorted = Vector{Vector{Float64}}(undef, J)
    x_sorted = Vector{Vector{Float64}}(undef, J)

    @threads for j in 1:J
        sub_obj_jth = @view sub_obj[:,j]; ascending_indices = sortperm(sub_obj_jth)
        sorted_indices[j] = ascending_indices; c_sorted[j] = sub_obj_jth[ascending_indices]; x_sorted[j] = x[ascending_indices]
    end
    return sorted_indices, c_sorted, x_sorted
end

function find_critical_item(c::Vector{Float64}, x::Vector{Float64})
    
    sum_x::Float64 = 0.0
    for (idx, val) in enumerate(x)
        sum_x += val
        if sum_x >= 1.0
            return idx
        end
    end
    throw(AlgorithmException("`k` cannot be `nothing` as sum(x) >= 2 is enforced. Check the models."))
end

function retrieve_dual_values(log, d, sorted_indices, c_sorted, critical_facility_indices)
    # Initialize
    log.dual_var[:δ] .= 0

    @threads for j in eachindex(d)
        c_sort_jth = c_sorted[j]; sort_idx_jth = sorted_indices[j]; k = critical_facility_indices[j]

        # Update σ
        log.dual_var[:σ][j] = c_sort_jth[k]

        # Update δ
        if k != 1
            c_critical = c_sort_jth[k]; δ_row = @view sort_idx_jth[1:k-1]; δ_ij = @view log.dual_var[:δ][δ_row, j]; c_before_critical = @view c_sort_jth[1:k-1]
            @. δ_ij = c_critical - c_before_critical
        end
    end
end 

function gen_figure(log, global_iter)

    x = 1:length(log.α_set)

    # Draw α
    plot(x, log.α_set,
    xlabel = "DD iteration",
    ylabel = "step size",
    legend = false)

    # savefig("Polyak_α.png")
    savefig("Polyak_α_yλorder_$(global_iter).png")
    # savefig("T1000x1000_Polyak_α_$(global_iter).png")
    
    # Draw LB
    # plot(x, log.obj_set,
    # xlabel = "Iteration",
    # ylabel = "LB",
    # ylims = (0, maximum(log.obj_set)),
    # legend = false)

    # # savefig("LB.png")
    # savefig("T1000x1000_Polyak_LB_$(global_iter).png")

    # # Draw res norm
    # plot(x, log.residual_norm_set,
    # xlabel = "Iteration",
    # ylabel = "Norm res",
    # legend = false)

    # savefig("Polyak_norm_res.png")
    # savefig("norm_res.png")

    # Draw λ_k norm
    plot(x, log.λ_k_diff_set,
    xlabel = "DD iteration",
    ylabel = "||opt_λ-λ||",
    legend = false)

    # savefig("Polyak_norm_λ_diff.png")
    savefig("Polyak_norm_λ_diff_yλorder_$(global_iter).png")
    # savefig("T1000x1000_Polyak_norm_λ_diff_$(global_iter).png")

    # plot(x, log.λ_k_norm_set,
    # xlabel = "Iteration",
    # ylabel = "Norm norm λ",
    # legend = false)

    # savefig("Polyak_norm_λ.png")

    # Draw trajectory of (||y^*-y^{approx}||, ||λ^*-λ^{approx}||)
    # plot(log.y_diff_set, log.λ_k_diff_set,
    # xlabel = "||y^*-y^k||",
    # ylabel = "||λ^*-λ^k||",
    # legend = false)

    # savefig("Polyak_trajectory(norm_y,norm_λ)_yλorder_$(global_iter).png")

end
