function build_cflp_core_point_x(data::CFLPData; optimizer = optimizer)
    model = Model(optimizer)
    I = data.n_facilities

    @variable(model, delta >= 0.0)
    @variable(model, x[1:I])
    @objective(model, Max, delta)

    @constraint(model, [i in 1:I], x[i] >= delta)
    @constraint(model, [i in 1:I], x[i] <= 1.0 - delta)
    @constraint(model, sum(data.capacities[i] * x[i] for i in 1:I) >= sum(data.demands))
    @constraint(model, sum(x) >= 1.0)

    optimize!(model)

    if termination_status(model) != OPTIMAL
        throw(UnexpectedModelStatusException("Unable to build CFLP core point: termination status $(termination_status(model))."))
    end

    return value.(x), objective_value(model)
end

build_cflp_ones_core_point_x(data::CFLPData) = (ones(data.n_facilities), 0.0)

function solve_cflp_recourse_value(data::CFLPData, x_value::Vector{Float64}; optimizer = optimizer)
    model = Model(optimizer)
    I, J = data.n_facilities, data.n_customers
    @variable(model, y[1:I, 1:J] >= 0)

    cost_demands = data.costs .* data.demands'
    @objective(model, Min, sum(cost_demands .* y))
    @constraint(model, demand[j in 1:J], sum(y[:, j]) == 1)
    @constraint(model, facility_open[i in 1:I, j in 1:J], y[i, j] <= x_value[i])
    @constraint(model, capacity[i in 1:I], sum(data.demands .* y[i, :]) <= data.capacities[i] * x_value[i])

    optimize!(model)
    if termination_status(model) != OPTIMAL
        return Inf
    end
    return objective_value(model)
end

function build_anchor_t(recourse_value::Float64, t_margin_rel::Float64, t_margin_abs::Float64)
    t_margin = max(t_margin_abs, t_margin_rel * max(1.0, recourse_value))
    return [recourse_value + t_margin], t_margin
end

function build_discorvernorm_support_points(
    data::CFLPData;
    t_margin_rel::Float64 = 0.0,
    t_margin_abs::Float64 = 0.0,
    optimizer = optimizer,
)
    ones_x, _ = build_cflp_ones_core_point_x(data)
    ones_recourse = solve_cflp_recourse_value(data, ones_x; optimizer = optimizer)
    isfinite(ones_recourse) || throw(UnexpectedModelStatusException("Unable to construct a finite recourse value for the ones support anchor."))
    ones_t, ones_margin = build_anchor_t(ones_recourse, t_margin_rel, t_margin_abs)

    zero_x = zeros(data.n_facilities)
    zero_t = [1.0]

    return (
        support_points_x = hcat(ones_x, zero_x),
        support_points_t = hcat(ones_t, zero_t),
        ones_recourse = ones_recourse,
        ones_margin = ones_margin,
        zero_t = zero_t[1],
    )
end

const DISCORNORM_FLCAP_INSTANCE_NAMES = [
    "cap61", "cap62", "cap63", "cap64",
    "cap71", "cap72", "cap73", "cap74",
    "cap91", "cap92", "cap93", "cap94",
    "cap101", "cap102", "cap103", "cap104",
    "cap121", "cap122", "cap123", "cap124",
    "cap131", "cap132", "cap133", "cap134",
    "capa1", "capa2", "capa3", "capa4",
    "capb1", "capb2", "capb3", "capb4",
    "capc1", "capc2", "capc3", "capc4",
]

function parse_instance_token(token::AbstractString)
    token = strip(token)
    isempty(token) && return Int[]
    if occursin(":", token)
        bounds = split(token, ":")
        length(bounds) == 2 || throw(ArgumentError("Unsupported instance range token `$(token)`."))
        lo = parse(Int, strip(bounds[1]))
        hi = parse(Int, strip(bounds[2]))
        lo <= hi || throw(ArgumentError("Invalid decreasing instance range `$(token)`."))
        return collect(lo:hi)
    end
    return [parse(Int, token)]
end

function discorvernorm_benchmark_instances()
    raw = get(ENV, "DISCORNORM_CFL_INSTANCES", "")
    if isempty(strip(raw))
        return setdiff(1:71, [67])
    end

    parsed = Int[]
    for token in split(raw, ",")
        append!(parsed, parse_instance_token(token))
    end
    unique!(parsed)
    sort!(parsed)
    return parsed
end

function discorvernorm_flcap_instances()
    raw = get(ENV, "DISCORNORM_FLCAP_INSTANCES", "")
    if isempty(strip(raw))
        return copy(DISCORNORM_FLCAP_INSTANCE_NAMES)
    end

    parsed = [strip(token) for token in split(raw, ",") if !isempty(strip(token))]
    invalid = setdiff(parsed, DISCORNORM_FLCAP_INSTANCE_NAMES)
    isempty(invalid) || throw(ArgumentError("Unsupported FLCAP instance names: $(join(invalid, ", "))."))
    return parsed
end
