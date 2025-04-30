
"""
Run DCGLP cutting-plane
"""

function solve_dcglp!(oracle::DisjunctiveOracle, x_value::Vector{Float64}, t_value::Vector{Float64}, zero_indices::Vector{Int64}, one_indices::Vector{Int64}; zero_tol = 1e-9, time_limit = time_limit, throw_typical_cuts_for_errors = true, include_disjuctive_cuts_to_hyperplanes = true)
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
                try 
                    optimize!(dcglp)
                catch e
                    if throw_typical_cuts_for_errors
                        @warn "Returning typical Benders cuts due to unexpected error encountered when optimizing dcglp master: $e"
                        return generate_cuts(typical_oracles[1], x_value, t_value; time_limit = get_sec_remaining(log.start_time, time_limit))
                    else
                        throw(UnexpectedModelStatusException("DCGLP master: unexpected error encountered when optimizing dcglp master: $e"))
                    end
                end
                if is_solved_and_feasible(dcglp; allow_local = false, dual = true)
                    for i=1:2
                        state.values[:ω_x][i] = value.(dcglp[:omega_x][i,:])
                        state.values[:ω_t][i] = value.(dcglp[:omega_t][i,:])
                        state.values[:ω_0][i] = value(dcglp[:omega_0][i])
                    end
                    state.values[:tau] = value(dcglp[:tau])
                    state.LB = state.values[:tau]
                elseif termination_status(dcglp) == ALMOST_INFEASIBLE
                    if throw_typical_cuts_for_errors
                        @warn "Returning typical Benders cuts due to unexpected dcglp master termination status: $(termination_status(dcglp)); the problem is infeasible or dcglp encountered numerical issue"
                        return generate_cuts(typical_oracles[1], x_value, t_value; time_limit = get_sec_remaining(log.start_time, time_limit))
                    else
                        throw(UnexpectedModelStatusException("DCGLP master: unexpected dcglp master termination status: $(termination_status(dcglp)); the problem is infeasible or dcglp encountered numerical issue"))
                    end
                elseif termination_status(dcglp) == TIME_LIMIT
                    throw(TimeLimitException("Time limit reached during dcglp solving"))
                else
                    throw(UnexpectedModelStatusException("DCGLP master: $(termination_status(dcglp))"))
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
                        state.is_in_L[i], hyperplanes_a, state.f_x[i] = generate_cuts(typical_oracles[i], clamp.(ω_x[i] / ω_0[i], 0.0, 1.0), ω_t[i] / ω_0[i], tol = zero_tol / ω_0[i], time_limit = get_sec_remaining(log.start_time, time_limit))
                        
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
                if typeof(oracle.oracle_param.norm) <: LpNorm
                    state.values[:sx] = value.(dcglp[:sx])
                    update_upper_bound_and_gap!(state, log, (t1, t2) -> LinearAlgebra.norm([state.values[:sx]; t1 .+ t2 .- t_value], oracle.oracle_param.norm.p))
                elseif typeof(oracle.oracle_param.norm) <: StandardNorm
                    # update_upper_bound_and_gap!(state, log, (t1, t2) -> LinearAlgebra.norm([state.values[:sx]; t1 .+ t2 .- t_value], oracle.oracle_param.norm.p))
                end
            end

            record_iteration!(log, state)
        end
        
        oracle.param.verbose && print_iteration_info(state, log, oracle.oracle_param.norm)

        check_lb_improvement!(state, log; zero_tol = zero_tol)

        is_terminated(state, log, oracle.param, time_limit) && break

        typeof(oracle.oracle_param.norm) <: LpNorm ? add_constraints(dcglp, :con_benders, [benders_cuts[1]; benders_cuts[2]]) : add_constraints(dcglp, :con_benders, [benders_cuts[1]; benders_cuts[2]]; lhs = dcglp[:tau]) 
    end

    if log.iterations[end].LB >= zero_tol
        if oracle.oracle_param.lift 
            gamma_x, gamma_t, gamma_0 = generate_lifted_disjunctive_cut(oracle.dcglp, oracle.oracle_param.norm, zero_indices, one_indices; strengthen = oracle.oracle_param.strengthened)
        else
            gamma_x, gamma_t, gamma_0 = generate_disjunctive_cut(oracle.dcglp; strengthen = oracle.oracle_param.strengthened)
        end

        h = Hyperplane(gamma_x, gamma_t, gamma_0)
        if include_disjuctive_cuts_to_hyperplanes
            push!(hyperplanes, h)
        end
        
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
            typeof(oracle.oracle_param.norm) <: LpNorm ? add_constraints(dcglp, :con_disjunctive, d_cuts)  : add_constraints(dcglp, :con_disjunctive, d_cuts; lhs = dcglp[:tau]) 
        end
        
        return false, hyperplanes, fill(Inf, length(t_value))
    else
        return generate_cuts(typical_oracles[1], x_value, t_value; time_limit = get_sec_remaining(log.start_time, time_limit))
    end
    # statistics_of_disjunctive_cuts(env)
end

function generate_disjunctive_cut(dcglp::Model; strengthen = false, zero_tol = 1e-7)
    gamma_x = dual.(dcglp[:conx])
    gamma_t = dual.(dcglp[:cont])
    gamma_0 = dual(dcglp[:con0])

    if strengthen 
        sigma = Dict(1 => dual(dcglp[:con_split_kappa]), 2 => dual(dcglp[:con_split_nu]))
        delta = Dict(1 => dual.(dcglp[:condelta][1,:]), 2 => dual.(dcglp[:condelta][2,:]))
        gamma_x = -strengthening!(-gamma_x, sigma, delta; zero_tol = zero_tol)
    end
    
    return gamma_x, gamma_t, gamma_0
end

function generate_lifted_disjunctive_cut(dcglp::Model, norm::LpNorm, zero_indices::Vector{Int64}, one_indices::Vector{Int64}; strengthen = false, zero_tol = 1e-7)
    gamma_x = dual.(dcglp[:conx])
    gamma_t = dual.(dcglp[:cont])
    gamma_0 = dual(dcglp[:con0])

    zeta_k = !isempty(zero_indices) ? dual.(dcglp[:con_zeta][1,:]) : Float64[]
    zeta_v = !isempty(zero_indices) ? dual.(dcglp[:con_zeta][2,:]) : Float64[] 
    xi_k = !isempty(one_indices) ? dual.(dcglp[:con_xi][1,:]) : Float64[] 
    xi_v = !isempty(one_indices) ? dual.(dcglp[:con_xi][2,:]) : Float64[] 

    # coefficients for lifted cut
    lifted_gamma_0 = gamma_0 - sum(max.(xi_k, xi_v))
    lifted_gamma_x = zeros(Float64, length(gamma_x))
    lifted_gamma_x .= -gamma_x

    lifted_gamma_x[zero_indices] = -gamma_x[zero_indices] .+ max.(zeta_k, zeta_v)
    lifted_gamma_x[one_indices] = -gamma_x[one_indices] .- max.(xi_k, xi_v)

    if strengthen
        sigma = Dict(1 => dual(dcglp[:con_split_kappa]), 2 => dual(dcglp[:con_split_nu]))
        delta_1 = dual.(dcglp[:condelta][1,:])
        delta_2 = dual.(dcglp[:condelta][2,:])
        delta_1[zero_indices] += (-zeta_k + max.(zeta_k, zeta_v))
        delta_2[zero_indices] += (-zeta_v + max.(zeta_k, zeta_v)) 
        delta = Dict(1 => delta_1, 2 => delta_2)
        lifted_gamma_x = strengthening!(lifted_gamma_x, sigma, delta; zero_tol = zero_tol)
    end

    # compute normalization value
    norm_value = compute_norm_value(lifted_gamma_x, gamma_t, norm)

    return (-lifted_gamma_x, gamma_t, lifted_gamma_0) ./ norm_value
end

function strengthening!(gamma_x, sigma, delta; zero_tol = 1e-7)
    @debug "dcglp strengthening - sigma values: [σ₁: $(sigma[1]), σ₂: $(sigma[2])]"
    @debug "dcglp strengthening - delta values: [δ₁: $(delta[1]), δ₂: $(delta[2])]"
    
    a₁ = gamma_x .- delta[1]
    a₂ = gamma_x .- delta[2]
    sigma_sum = sigma[1] + sigma[2]
    if sigma_sum >= zero_tol
        m = (a₁ .- a₂) / sigma_sum
        m_lb = floor.(m)
        m_ub = ceil.(m)
        gamma_x = min.(a₁-sigma[1]*m_lb, a₂+sigma[2]*m_ub)
    end
    return gamma_x
end

function compute_norm_value(gamma_x, gamma_t, norm::AbstractNorm)
    # compute normalization value
    if norm.p == 1.0
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), Inf)
    elseif norm.p == 2.0
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), 2.0)
    elseif norm.p == Inf
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), 1.0)
    else
        throw(UndefError("Unsupported norm type: $(typeof(norm))"))
    end
    return max(1.0, norm_value)
end