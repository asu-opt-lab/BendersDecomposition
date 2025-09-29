export update_λ, update_y, retrieve_dual_values, gen_figure

"""
utility functions for Dual Decomposition
"""
function update_λ(λ_k, y_k, u, d, x, α)
    for i in eachindex(λ_k)
        λ_k[i] = max(0, λ_k[i] + α*(d'*y_k[i,:] - u[i]*x[i]))
    end    
    return λ_k
end

function update_y(y_k, λ_k, c, d, x)
    J = collect(1:length(d))
    critical_facility_indices = Vector{Int}(undef, length(J))
    sorted_indices, c_sorted, x_sorted = sort_costs_indices(J, c + λ_k .* d', x) # sorted_indice(x index): sort x with cheapest order, x_sorted (x value sorted)

    for j in J
        k = find_critical_item(c_sorted[j], x_sorted[j]) # k: index of critical item, sorted_indices[j][k] is a critical item
        critical_facility_indices[j] = k
        jth_sortin = sorted_indices[j]

        y_col = y_k[:,j] # choose jth column and modify row values as follows

        # k != 1 && (y_col[jth_sortin[1:k-1]] .= x_sorted[j][k-1])
        y_col[jth_sortin[1:k-1]] .= k != 1 ? x_sorted[j][k-1] : y_col[jth_sortin[1:k-1]] # LHS: indices < k, RHS: corresponidng x values
        y_col[jth_sortin[k]] = k != 1 ? 1 - sum(x_sorted[j][k-1]) : x_sorted[j][k] # if index of critical item is not 1, some y > 0 before k. Assign residual value to corresponding y.
        y_col[jth_sortin[k+1:end]] .= 0 # facility index > index of critical item should be all zero

        y_k[:, j] .= y_col
    end
    return y_k, sorted_indices, c_sorted, x_sorted, critical_facility_indices
end

function sort_costs_indices(J, sub_obj, x)
    sorted_indices, c_sorted, x_sorted = [], [], []
    for j in J
        push!(sorted_indices, sortperm(sub_obj[:,j]))
        push!(c_sorted, sub_obj[:,j][sorted_indices[j]]); push!(x_sorted, x[sorted_indices[j]])
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

    for j in eachindex(d)
        c_sort_jth = c_sorted[j]; sort_idx_jth = sorted_indices[j]; k = critical_facility_indices[j]

        # Update σ
        log.dual_var[:σ][j] = c_sort_jth[k]

        # Update δ
        k != 1 && (log.dual_var[:δ][sort_idx_jth[1:k-1],j] .= c_sort_jth[k] .- c_sort_jth[1:k-1])
    end
end 

function gen_figure(log, global_iter)

    x = 1:length(log.α_set)

    # Draw α
    plot(x, log.α_set,
    xlabel = "Iteration",
    ylabel = "step size",
    legend = false)

    # savefig("Polyak_α.png")
    # savefig("α.png")
    savefig("T1000x1000_Polyak_α_$(global_iter).png")
    

    # Draw LB
    plot(x, log.obj_set,
    xlabel = "Iteration",
    ylabel = "LB",
    ylims = (0, maximum(log.obj_set)),
    legend = false)

    # savefig("LB.png")
    savefig("T1000x1000_Polyak_LB_$(global_iter).png")

    # # Draw obj without penalty
    # plot(x, log.approx_obj_set,
    # xlabel = "Iteration",
    # ylabel = "sum(c.*y)",
    # ylims = (0, maximum(log.approx_obj_set)),
    # legend = false)

    # savefig("Polyak_sum(c.*y).png")
    # savefig("obj_sum(c.*y).png")
    # savefig("T1000x1000_Polyak_sum(c.*y).png")

    # # Draw res norm
    # plot(x, log.residual_norm_set,
    # xlabel = "Iteration",
    # ylabel = "Norm res",
    # legend = false)

    # savefig("Polyak_norm_res.png")
    # savefig("norm_res.png")

    # Draw λ_k norm
    plot(x, log.λ_k_diff_set,
    xlabel = "Iteration",
    ylabel = "||opt_λ-λ||",
    legend = false)

    # savefig("Polyak_norm_λ_diff.png")
    savefig("T1000x1000_Polyak_norm_λ_diff_$(global_iter).png")

    # plot(x, log.λ_k_norm_set,
    # xlabel = "Iteration",
    # ylabel = "Norm norm λ",
    # legend = false)

    # savefig("Polyak_norm_λ.png")

end
