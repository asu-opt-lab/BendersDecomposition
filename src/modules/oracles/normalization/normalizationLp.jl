# -----------------------------------------------------------------------------
# Lp-distance normalization
# -----------------------------------------------------------------------------

function add_normalization_constraint!(
    normalization::LpDistanceNormalization,
    dcglp::Model,
    tau::VariableRef,
    sx::AbstractVector{VariableRef},
    st::AbstractVector{VariableRef},
)

    var_vec = [tau; sx; st]
    norm = normalization.norm_p
    norm isa LpNorm || throw(UndefError("Unsupported norm type: $(typeof(norm))"))
    if norm.p == 1.0
        @constraint(dcglp, concone, var_vec in MOI.NormOneCone(length(var_vec)))
    elseif norm.p == 2.0
        @constraint(dcglp, concone, var_vec in MOI.SecondOrderCone(length(var_vec)))
    elseif norm.p == Inf
        @constraint(dcglp, concone, var_vec in MOI.NormInfinityCone(length(var_vec)))
    else
        throw(UndefError("Unsupported LpNorm: p=$(norm.p)"))
    end
end

function update_dcglp_upper_bound_and_gap!(
    normalization::LpDistanceNormalization,
    state::DcglpState,
    log::DcglpLog,
    reference_t::Vector{Float64},
    t_value::Vector{Float64},
)
    fill_dcglp_omega_t_estimates!(state, t_value)
    all(f_i -> !any(isnan, f_i), state.f_x) || return nothing
    update_upper_bound_and_gap!(
        state,
        log,
        (t1, t2) -> LinearAlgebra.norm([state.values[:sx]; t1 .+ t2 .- reference_t], normalization.norm_p.p),
    )
end

function compute_norm_value(gamma_x, gamma_t, norm::LpNorm)
    if norm.p == 1.0
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), Inf)
    elseif norm.p == 2.0
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), 2.0)
    elseif norm.p == Inf
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), 1.0)
    else
        throw(UndefError("Unsupported LpNorm: p=$(norm.p)"))
    end
    return max(1.0, norm_value)
end

function build_dcglp_disjunctive_cut(
    normalization::LpDistanceNormalization,
    dcglp::Model,
    common::SplitOracleParam,
    ::Float64,
    ::Vector{Float64},
    ::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int},
)
    gamma_x = dual.(dcglp[:conx])
    gamma_t = dual.(dcglp[:cont])
    gamma_0 = dual(dcglp[:con0])

    gamma_x, gamma_0 = apply_lift_or_strengthen(
        dcglp, gamma_x, zero_indices, one_indices;
        lift = common.lift, strengthen = common.strengthened,
        zero_tol = common.zero_tol, gamma_0 = gamma_0,
    )

    if common.lift && (!isempty(zero_indices) || !isempty(one_indices))
        norm_value = compute_norm_value(gamma_x, gamma_t, normalization.norm_p)
        return Hyperplane(gamma_x ./ norm_value, gamma_t ./ norm_value, gamma_0 / norm_value)
    end

    return Hyperplane(gamma_x, gamma_t, gamma_0)
end
