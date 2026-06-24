"""
    AbstractDisjunctiveNormalization

Normalization object that captures the per-variant behavior of a [`SplitOracle`](@ref).
All split-oracle normalization variants in this module are stored on `SplitOracleParam`.
The normalization owns the variant-specific
configuration (`norm`, `core_point_*`, …) and the dispatch surface listed below.

# Normalization interface
- `build_dcglp` — construct the DCGLP for a normalization. Normalizations with
  distinct linking constraints provide their own method.
- `add_normalization_constraint!` — connect the shared `tau`, `sx`, and `st`
  variables according to the normalization.
- `update_dcglp_for_candidate!` — per-call setup for the current
  candidate point (RHS, direction).
- `update_dcglp_reference_t!` — normalizations may override to adjust the
  reference epigraph used by the UB recomputation.
- `update_dcglp_upper_bound_and_gap!` — normalization-specific UB/gap rule.
- `disjunctive_cut_normalization_value` — normalization-specific scaling in the common
  disjunctive cut extraction.
- `print_dcglp_iteration_info` — normalization-specific iteration log.
"""
abstract type AbstractDisjunctiveNormalization end

# Default interface implementations (overridden per normalization where needed).

function update_dcglp_for_candidate!(::AbstractDisjunctiveNormalization, oracle::AbstractSplitOracle, x_value::Vector{Float64}, t_value::Vector{Float64})
    set_normalized_rhs.(oracle.dcglp[:conx], x_value)
    set_normalized_rhs.(oracle.dcglp[:cont], t_value)
    return nothing
end

function prepare_disjunctive_normalization!(::AbstractDisjunctiveNormalization, ::AbstractMaster)
    return nothing
end

function add_normalization_constraint!(
    normalization::AbstractDisjunctiveNormalization,
    dcglp::Model,
    tau::VariableRef,
    sx::AbstractVector{VariableRef},
    st::AbstractVector{VariableRef},
)
    throw(UndefError("update add_normalization_constraint! for $(typeof(normalization))"))
end

function update_dcglp_upper_bound_and_gap!(
    normalization::AbstractDisjunctiveNormalization,
    state::DcglpState,
    log::DcglpLog,
    reference_t::Vector{Float64},
    t_value::Vector{Float64},
)
    throw(UndefError("update update_dcglp_upper_bound_and_gap! for $(typeof(normalization))"))
end

function adjust_dcglp_reference_t_to_fx!(
    oracle::AbstractSplitOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    start_time::Float64,
    time_limit::Float64,
)
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

function update_dcglp_reference_t!(
    ::AbstractDisjunctiveNormalization,
    oracle::AbstractSplitOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    start_time::Float64,
    time_limit::Float64,
)
    oracle.param.adjust_t_to_fx || return copy(t_value)
    return adjust_dcglp_reference_t_to_fx!(oracle, x_value, t_value, start_time, time_limit)
end
