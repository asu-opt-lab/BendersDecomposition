using BendersX
using JuMP
using CPLEX
using Printf
using Random
using MathOptInterface
using LinearAlgebra

isdefined(Main, :PolarDCGLP) || include(normpath(joinpath(@__DIR__, "..", "src", "PolarDCGLP.jl")))
using .PolarDCGLP

include(normpath(joinpath(@__DIR__, "solver_defaults.jl")))
include(normpath(joinpath(@__DIR__, "script_utils.jl")))

const MOI = MathOptInterface

function parse_split_rule(name::String)
    normalized = lowercase(strip(name))
    if normalized in ("most_fractional", "most", "mf")
        return MostFractional()
    elseif normalized in ("largest_fractional", "largest", "lf")
        return LargestFractional()
    elseif normalized in ("random_fractional", "random", "rf")
        return RandomFractional()
    end
    throw(ArgumentError("Unsupported split rule `$(name)`. Use one of: most_fractional, largest_fractional, random_fractional."))
end

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

function build_cflp_ones_core_point_x(data::CFLPData)
    return ones(data.n_facilities), 0.0
end

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

function build_directional_core_point(
    data::CFLPData;
    x_mode::String,
    t_margin_rel::Float64,
    t_margin_abs::Float64,
    optimizer = optimizer,
)
    normalized_x_mode = lowercase(strip(x_mode))
    if normalized_x_mode in ("ones", "all_ones", "one")
        x_core, delta = build_cflp_ones_core_point_x(data)
    elseif normalized_x_mode in ("centered", "interior", "max_delta")
        x_core, delta = build_cflp_core_point_x(data; optimizer = optimizer)
    else
        throw(ArgumentError("Unsupported core_x_mode `$(x_mode)`. Use one of: ones, centered."))
    end

    recourse_value = solve_cflp_recourse_value(data, x_core; optimizer = optimizer)

    if !isfinite(recourse_value)
        @warn "Auto-generated CFLP core point was not subproblem-feasible; falling back to x = 1."
        x_core .= 1.0
        recourse_value = solve_cflp_recourse_value(data, x_core; optimizer = optimizer)
    end

    isfinite(recourse_value) || throw(UnexpectedModelStatusException("Unable to construct a finite core-point recourse value for CFLP."))

    t_margin = max(t_margin_abs, t_margin_rel * max(1.0, recourse_value))
    return x_core, [recourse_value + t_margin], delta, recourse_value
end

options, _ = parse_script_args(ARGS)

instance = get_string_option(options, "instance", "capc3")
seed = get_int_option(options, "seed", 1)
time_limit = get_float_option(options, "time_limit", 1200.0)
dcglp_time_limit = get_float_option(options, "dcglp_time_limit", 1000.0)
dcglp_iter_limit = get_int_option(options, "dcglp_iter_limit", 100)
dcglp_halt_limit = get_int_option(options, "dcglp_halt_limit", 3)
frequency = get_int_option(options, "frequency", 2500)
reuse_dcglp = get_bool_option(options, "reuse_dcglp", false)
strengthened = get_bool_option(options, "strengthened", false)
build_only = get_bool_option(options, "build_only", false)
split_rule_name = get_string_option(options, "split_rule", "most_fractional")
core_x_mode = get_string_option(options, "core_x_mode", "ones")
core_t_margin_rel = get_float_option(options, "core_t_margin_rel", 0.05)
core_t_margin_abs = get_float_option(options, "core_t_margin_abs", 1.0)

Random.seed!(seed)

@info "DirectionalPolarDCGLP CFL script" instance = instance seed = seed time_limit = time_limit frequency = frequency reuse_dcglp = reuse_dcglp strengthened = strengthened build_only = build_only split_rule = split_rule_name core_x_mode = core_x_mode

data = read_flcap_data(instance)
split_rule = parse_split_rule(split_rule_name)
core_x, core_t, core_delta, core_recourse = build_directional_core_point(
    data;
    x_mode = core_x_mode,
    t_margin_rel = core_t_margin_rel,
    t_margin_abs = core_t_margin_abs,
    optimizer = optimizer,
)

@info "Directional core point constructed" instance = instance core_x_mode = core_x_mode core_delta = core_delta core_x_min = minimum(core_x) core_x_max = maximum(core_x) capacity_ratio = dot(data.capacities, core_x) / sum(data.demands) core_recourse = core_recourse core_t = core_t[1]

benders_param = BendersBnBParam(
    time_limit = time_limit,
    gap_tolerance = 1e-6,
    verbose = true,
)

dcglp_optimizer = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_Threads" => 7,
    "CPX_PARAM_EPRHS" => 1e-9,
    "CPX_PARAM_NUMERICALEMPHASIS" => 1,
    "CPX_PARAM_EPOPT" => 1e-9,
    MOI.Silent() => true,
)
dcglp_param = DcglpParam(
    dcglp_optimizer;
    time_limit = dcglp_time_limit,
    gap_tolerance = 1e-3,
    halt_limit = dcglp_halt_limit,
    iter_limit = dcglp_iter_limit,
    verbose = true,
)

oracle_param = DirectionalPolarDCGLPParam(
    dcglp_param,
    core_x,
    core_t;
    split_index_selection_rule = split_rule,
    disjunctive_cut_append_rule = AllDisjunctiveCuts(),
    add_benders_cuts_to_master = true,
    fraction_of_benders_cuts_to_master = 0.5,
    reuse_dcglp = reuse_dcglp,
    strengthened = strengthened,
    zero_tol = 1e-9,
)

master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)

typical_oracles = [
    CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
    CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
]
disjunctive_oracle = DirectionalPolarDCGLPOracle(master, typical_oracles, oracle_param)

lazy_oracle = CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer)
root_preprocessing = RootNodePreprocessing(
    lazy_oracle,
    BendersSeqInOut,
    BendersSeqInOutParam(
        time_limit = min(100.0, time_limit),
        gap_tolerance = 1e-6,
        stabilizing_x = ones(data.n_facilities),
        α = 0.9,
        λ = 0.1,
        verbose = true,
    ),
)
lazy_callback = LazyCallback(lazy_oracle)
user_callback = UserCallback(disjunctive_oracle; params = UserCallbackParam(frequency = frequency))

env = BendersBnB(
    master,
    root_preprocessing,
    lazy_callback,
    user_callback;
    param = benders_param,
)

if build_only
    @info "DirectionalPolarDCGLP CFL script build completed without solve." instance = instance
else
    solve!(env)
    obj_value = try
        env.obj_value
    catch
        NaN
    end
    @info "DirectionalPolarDCGLP CFL script finished" instance = instance termination_status = env.termination_status objective_value = obj_value
end
