
"""
Run DCGLP cutting-plane
"""
# lifting
function solve_dcglp!(oracle::DisjunctiveOracle, x_value::Vector{Float64}, t_value::Vector{Float64}, zero_indices::Vector{Int64}, one_indices::Vector{Int64}; zero_tol = 1e-9, time_limit = time_limit)
    log = DcglpLog()
    
    dcglp = oracle.dcglp
    typical_oracles = oracle.typical_oracles # 1 for k; 2 for nu

    # cuts for master
    hyperplanes = Vector{Hyperplane}()

    while true
        state = DcglpState()
        state.total_time = @elapsed begin
            state.master_time = @elapsed begin
                set_time_limit_sec(dcglp, get_sec_remaining(log.start_time, time_limit))
                optimize!(dcglp)
                if is_solved_and_feasible(dcglp; allow_local = false, dual = true)
                    for i=1:2
                        state.values[:ω_x][i] = value.(dcglp[:omega_x][i,:])
                        state.values[:ω_t][i] = value.(dcglp[:omega_t][i,:])
                        state.values[:ω_0][i] = value(dcglp[:omega_0][i])
                    end
                    state.values[:tau] = value(dcglp[:tau])
                    state.values[:sx] = value.(dcglp[:sx])
                    state.LB = state.values[:tau]
                elseif termination_status(dcglp) == ALMOST_INFEASIBLE
                    @warn "dcglp master termination status: $(termination_status(dcglp)); the problem is infeasible or dcglp encountered numerical issue, yielding the typical Benders cut"
                    return generate_cuts(typical_oracles[1], x_value, t_value; time_limit = get_sec_remaining(log.start_time, time_limit))
                elseif termination_status(dcglp) == TIME_LIMIT
                    throw(TimeLimitException("Time limit reached during dcglp solving"))
                else
                    throw(UnexpectedModelStatusException("dcglp master: $(termination_status(dcglp))"))
                    # if infeasible, then the problem is infeasible
                end
            end
            
            benders_cuts = Dict(1 => Vector{AffExpr}(), 2 => Vector{AffExpr}())
            ω_x = state.values[:ω_x]
            ω_t = state.values[:ω_t]
            ω_0 = state.values[:ω_0]
            # my_lock = Threads.ReentrantLock()
            # Threads.@threads for i in 1:2
            for i in 1:2
                # Threads.lock(my_lock) do
                state.oracle_times[i] = @elapsed begin
                    if ω_0[i] >= zero_tol
                        state.is_in_L[i], hyperplanes_a, state.f_x[i] = generate_cuts(typical_oracles[i], ω_x[i] / ω_0[i], ω_t[i] / ω_0[i], time_limit = get_sec_remaining(log.start_time, time_limit))
                        # adjust the tolerance with respect to dcglp: (sum(state.sub_obj_vals[i]) - sum(t_value)) * omega_value[:z][i] < zero_tol
                        if !state.is_in_L[i]
                            for k = 1:2 # add to both kappa and nu systems
                                append!(benders_cuts[i], hyperplanes_to_expression(dcglp, hyperplanes_a, dcglp[:omega_x][k,:], dcglp[:omega_t][k,:], dcglp[:omega_0][k]))
                            end
                            if oracle.oracle_param.add_benders_cuts_to_master
                                append!(hyperplanes, select_top_fraction(hyperplanes_a, h -> evaluate_violation(h, x_value, t_value), oracle.oracle_param.fraction_of_benders_cuts_to_master))
                            end
                        end
                    else
                        state.is_in_L[i] = true
                        state.f_x[i] = zeros(length(t_value))
                    end
                end
                # end
            end
        
            if state.f_x[1] !== NaN && state.f_x[2] !== NaN
                update_upper_bound_and_gap!(state, log, (t1, t2) -> LinearAlgebra.norm([state.values[:sx]; t1 .+ t2 .- t_value], oracle.oracle_param.norm.p))
            end

            record_iteration!(log, state)
            oracle.param.verbose && print_iteration_info(state, log)
        end

        check_lb_improvement!(state, log; zero_tol = zero_tol)

        is_terminated(state, log, oracle.param, time_limit) && break

        add_constraints(dcglp, :con_benders, [benders_cuts[1]; benders_cuts[2]]) 
    end

    if log.iterations[end].LB >= zero_tol
        # lifting: we can shorten it by merging the functions and receive argument (lift), but the merged function will become longer
        if oracle.oracle_param.lift == true
            gamma_x, gamma_t, gamma_0 = oracle.oracle_param.strengthened ? generate_strengthened_lifted_disjunctive_cuts(oracle.dcglp, zero_indices, one_indices) : generate_lifted_disjunctive_cut(oracle.dcglp, zero_indices, one_indices)
        else
            gamma_x, gamma_t, gamma_0 = oracle.oracle_param.strengthened ? generate_strengthened_disjunctive_cuts(oracle.dcglp) : generate_disjunctive_cut(oracle.dcglp)
        end

        h = Hyperplane(gamma_x, gamma_t, gamma_0)
        push!(hyperplanes, h)
        
        if typeof(oracle.oracle_param.split_index_selection_rule) <: SimpleSplit
            index = get_split_index(oracle)
            push!(oracle.disjunctiveCutsByIndex[index], h)
        end
        push!(oracle.disjunctiveCuts, h)
        
        if oracle.oracle_param.disjunctive_cut_append_rule == AllDisjunctiveCuts()
            d_cuts = Vector{AffExpr}()
            for k = 1:2 # add to both kappa and nu systems
                append!(d_cuts, hyperplanes_to_expression(dcglp, [h], dcglp[:omega_x][k,:], dcglp[:omega_t][k,:], dcglp[:omega_0][k]))
            end
            add_constraints(dcglp, :con_disjunctive, d_cuts) 
        end
        
        return false, hyperplanes, fill(Inf, length(t_value))
    else
        return generate_cuts(typical_oracles[1], x_value, t_value; time_limit = get_sec_remaining(log.start_time, time_limit))
    end
    # statistics_of_disjunctive_cuts(env)
end

function generate_disjunctive_cut(dcglp::Model)
    gamma_x = dual.(dcglp[:conx])
    gamma_t = dual.(dcglp[:cont])
    gamma_0 = dual(dcglp[:con0])
    
    return gamma_x, gamma_t, gamma_0
end

function generate_strengthened_disjunctive_cuts(dcglp::Model; zero_tol = 1e-5)
    
    σ₁ = dual(dcglp[:con_split_kappa])
    σ₂ = dual(dcglp[:con_split_nu])
    gamma_x = dual.(dcglp[:conx])
    gamma_t = dual.(dcglp[:cont])
    gamma_0 = dual(dcglp[:con0])

    # println("DCGLP Sigma Values: [σ₁: $σ₁, σ₂: $σ₂]")

    a₁ = -gamma_x .- dual.(dcglp[:condelta][1])
    a₂ = -gamma_x .- dual.(dcglp[:condelta][2])
    σ_sum = σ₂ + σ₁
    if σ_sum >= zero_tol
        m = (a₁ .- a₂) / σ_sum
        m_lb = floor.(m)
        m_ub = ceil.(m)
        gamma_x = -min.(a₁-σ₁*m_lb, a₂+σ₂*m_ub)
    end
    
    return gamma_x, gamma_t, gamma_0
end

function generate_lifted_disjunctive_cut(oracle::DisjunctiveOracle, zero_indices::Vector{Float64}, one_indices::Vector{Float64})
    # we can remove gamma_x,t,0 if we use generate_disjunctive_cut default outside this function.
    gamma_x, gamma_t, gamma_0 = generate_disjunctive_cut(oracle.dcglp)
    @debug gamma_x
    @debug gamma_t
    @debug gamma_0
    error("bp")
    zeta_k = dual.(dcglp.model[:con_zeta][1]) # zero indices
    zeta_v = dual.(dcglp.model[:con_zeta][2]) # zero indices
    xi_k = dual.(dcglp.model[:con_xi][1]) # one indices
    xi_v = dual.(dcglp.model[:con_xi][2]) # one indices

    lifted_gamma_0 = gamma_0 - sum(max.(xi_k, xi_v))
    lifted_gamma_x = zeros(Float64, length(gamma_x))
    lifted_gamma_x .= -gamma_x

    for (i, idx) in enumerate(zeros_indices)
        lifted_gamma_x[idx] = -gamma_x[idx] + max(zeta_k[i], zeta_v[i])
    end
    for (i, idx) in enumerate(ones_indices)
        lifted_gamma_x[idx] = -gamma_x[idx] - max(xi_k[i],xi_v[i])
    end
    # lift_gamma_x[zeros_indices] = -gamma_x[zeros_indices] .+ max.(zeta_k, zeta_v)
    # lifted_gamma_x[ones_indices] = -gamma_x[ones_indices] .- max.(xi_k, xi_v)

    if norm_type == L1Norm()
        _norm_value = norm(vcat(lifted_gamma_x, gamma_t), 1)
    elseif norm_type == L2Norm()
        _norm_value = norm(vcat(lifted_gamma_x, gamma_t), 2)
    else
        _norm_value = norm(vcat(lifted_gamma_x, gamma_t), Inf)
    end
    # norm_value = norm(vcat(lifted_gamma_x, gamma_t), norm_type == L1Norm() ? 1 : norm_type == L2Norm() ? 2 : norm_type == LInfNorm() ? Inf : error("Unsupported norm_type: $norm_type"))
    norm_value = max(1.0, _norm_value)

    lifted_gamma_x = -(lifted_gamma_x ./ norm_value)
    gamma_t = gamma_t ./ norm_value
    lifted_gamma_0 = lifted_gamma_0 ./ norm_value


    return lifted_gamma_x, gamma_t, lifted_gamma_0
end


