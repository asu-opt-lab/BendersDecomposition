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
    norm = normalization.norm
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

function update_dcglp_reference_t!(
    ::LpDistanceNormalization,
    oracle::SplitOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    start_time::Float64,
    time_limit::Float64,
)
    reference_t = copy(t_value)
    oracle.param.normalization.adjust_t_to_fx || return reference_t

    dcglp = oracle.dcglp
    delete_registered_constraints!(dcglp, :initial_L)

    _, initial_hyperplanes, f_x = generate_cuts(
        oracle.typical_oracles[1],
        x_value,
        t_value;
        time_limit = get_sec_remaining(start_time, time_limit),
    )
    any(isnan, f_x) &&
        throw(AlgorithmException("solve_dcglp!: `t_value` cannot be adjusted to `f(x)` since $(typeof(oracle.typical_oracles[1])) does not compute `f(x)`."))

    initial_benders_cuts = AffExpr[]
    for k in 1:2
        append!(
            initial_benders_cuts,
            hyperplanes_to_expression(
                dcglp,
                initial_hyperplanes,
                dcglp[:omega_x][k, :],
                dcglp[:omega_t][k, :],
                dcglp[:omega_0][k],
            ),
        )
    end
    dcglp[:initial_L] = @constraint(dcglp, 0 .>= initial_benders_cuts)
    set_normalized_rhs.(dcglp[:cont], f_x)
    return copy(f_x)
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
        (t1, t2) -> LinearAlgebra.norm([state.values[:sx]; t1 .+ t2 .- reference_t], normalization.norm.p),
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
        norm_value = compute_norm_value(gamma_x, gamma_t, normalization.norm)
        return Hyperplane(gamma_x ./ norm_value, gamma_t ./ norm_value, gamma_0 / norm_value)
    end

    return Hyperplane(gamma_x, gamma_t, gamma_0)
end
