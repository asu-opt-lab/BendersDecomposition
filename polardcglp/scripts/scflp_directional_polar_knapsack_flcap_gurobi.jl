using BendersX
using JuMP
using Gurobi
using Printf
using Random
using MathOptInterface
using LinearAlgebra
using Statistics

const MOI = MathOptInterface

isdefined(Main, :PolarDCGLP) || include(normpath(joinpath(@__DIR__, "..", "src", "PolarDCGLP.jl")))
using .PolarDCGLP

include(normpath(joinpath(@__DIR__, "script_utils.jl")))

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

function build_scflp_ones_core_point_x(data::SCFLPData)
    return ones(data.n_facilities), 0.0
end

function build_scflp_centered_core_point_x(data::SCFLPData; optimizer)
    model = Model(optimizer)
    I = data.n_facilities

    @variable(model, delta >= 0.0)
    @variable(model, x[1:I])
    @objective(model, Max, delta)

    @constraint(model, [i in 1:I], x[i] >= delta)
    @constraint(model, [i in 1:I], x[i] <= 1.0 - delta)

    max_demand = maximum(sum(d) for d in data.demands)
    @constraint(model, sum(data.capacities[i] * x[i] for i in 1:I) >= max_demand)
    @constraint(model, sum(x) >= 1.0)

    optimize!(model)

    if termination_status(model) != OPTIMAL
        throw(UnexpectedModelStatusException("Unable to build SCFLP centered core point: termination status $(termination_status(model))."))
    end

    return value.(x), objective_value(model)
end

function solve_scflp_recourse_value(data::SCFLPData, x_value::Vector{Float64}, scen_idx::Int; optimizer)
    model = Model(optimizer)
    I, J = data.n_facilities, data.n_customers
    d_k = data.demands[scen_idx]

    @variable(model, y[1:I, 1:J] >= 0)

    cost_demands = data.costs .* d_k'
    @objective(model, Min, sum(cost_demands .* y))
    @constraint(model, demand[j in 1:J], sum(y[:, j]) == 1)
    @constraint(model, facility_open[i in 1:I, j in 1:J], y[i, j] <= x_value[i])
    @constraint(model, capacity[i in 1:I], sum(d_k .* y[i, :]) <= data.capacities[i] * x_value[i])

    optimize!(model)
    if termination_status(model) != OPTIMAL
        return Inf
    end
    return objective_value(model)
end

function build_scflp_directional_core_point(
    data::SCFLPData;
    x_mode::String,
    t_margin_rel::Float64,
    t_margin_abs::Float64,
    optimizer,
)
    normalized_x_mode = lowercase(strip(x_mode))
    if normalized_x_mode in ("ones", "all_ones", "one")
        x_core, delta = build_scflp_ones_core_point_x(data)
    elseif normalized_x_mode in ("centered", "interior", "max_delta")
        x_core, delta = build_scflp_centered_core_point_x(data; optimizer = optimizer)
    else
        throw(ArgumentError("Unsupported core_x_mode `$(x_mode)`. Use one of: ones, centered."))
    end

    recourse_per_scen = [solve_scflp_recourse_value(data, x_core, k; optimizer = optimizer) for k in 1:data.n_scenarios]

    if any(!isfinite, recourse_per_scen)
        @warn "Auto-generated SCFLP core point was not subproblem-feasible for some scenarios; falling back to x = 1."
        x_core = ones(data.n_facilities)
        recourse_per_scen = [solve_scflp_recourse_value(data, x_core, k; optimizer = optimizer) for k in 1:data.n_scenarios]
    end

    all(isfinite, recourse_per_scen) || throw(UnexpectedModelStatusException("Unable to construct finite core-point recourse values for SCFLP."))

    core_t = Vector{Float64}(undef, data.n_scenarios)
    for k in 1:data.n_scenarios
        margin = max(t_margin_abs, t_margin_rel * max(1.0, recourse_per_scen[k]))
        core_t[k] = recourse_per_scen[k] + margin
    end

    return x_core, core_t, delta, mean(recourse_per_scen)
end

options, _ = parse_script_args(ARGS)

instance = get_string_option(options, "instance", "cap101-s256")
seed = get_int_option(options, "seed", 1)
time_limit = get_float_option(options, "time_limit", 14400.0)
dcglp_time_limit = get_float_option(options, "dcglp_time_limit", 1000.0)
dcglp_iter_limit = get_int_option(options, "dcglp_iter_limit", 250)
dcglp_halt_limit = get_int_option(options, "dcglp_halt_limit", 3)
frequency = get_int_option(options, "frequency", 250)
threads = get_int_option(options, "threads", 7)
reuse_dcglp = get_bool_option(options, "reuse_dcglp", false)
strengthened = get_bool_option(options, "strengthened", true)
build_only = get_bool_option(options, "build_only", false)
split_rule_name = get_string_option(options, "split_rule", "most_fractional")
core_x_mode = get_string_option(options, "core_x_mode", "ones")
core_t_margin_rel = get_float_option(options, "core_t_margin_rel", 0.05)
core_t_margin_abs = get_float_option(options, "core_t_margin_abs", 1.0)

Random.seed!(seed)

@info "DirectionalPolarDCGLP SCFLP script (Gurobi, FLCAP)" instance = instance seed = seed time_limit = time_limit frequency = frequency threads = threads reuse_dcglp = reuse_dcglp strengthened = strengthened build_only = build_only split_rule = split_rule_name core_x_mode = core_x_mode

# Gurobi-based mip_optimizer (replaces solver_defaults.jl which uses CPLEX)
mip_optimizer = optimizer_with_attributes(
    Gurobi.Optimizer,
    "Threads" => threads,
    "IntFeasTol" => 1e-9,
    "FeasibilityTol" => 1e-9,
    "MIPGap" => 1e-6,
    "OptimalityTol" => 1e-9,
    "NumericFocus" => 1,
    MOI.Silent() => true,
)

optimizer = optimizer_with_attributes(
    Gurobi.Optimizer,
    "Threads" => threads,
    "FeasibilityTol" => 1e-9,
    "OptimalityTol" => 1e-9,
    "NumericFocus" => 1,
    MOI.Silent() => true,
)

dcglp_optimizer = optimizer_with_attributes(
    Gurobi.Optimizer,
    "Threads" => threads,
    "FeasibilityTol" => 1e-9,
    "NumericFocus" => 1,
    "OptimalityTol" => 1e-9,
    MOI.Silent() => true,
)

const SCFLP_DIR = normpath(joinpath(@__DIR__, "..", "data", "SCFLP"))

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
data = read_stochastic_capacited_facility_location_problem(instance; filepath = SCFLP_DIR)

split_rule = parse_split_rule(split_rule_name)

core_x, core_t, core_delta, mean_recourse = build_scflp_directional_core_point(
    data;
    x_mode = core_x_mode,
    t_margin_rel = core_t_margin_rel,
    t_margin_abs = core_t_margin_abs,
    optimizer = optimizer,
)

@info "Directional core point constructed" instance = instance core_x_mode = core_x_mode core_delta = core_delta core_x_min = minimum(core_x) core_x_max = maximum(core_x) capacity_ratio = dot(data.capacities, core_x) / maximum(sum(d) for d in data.demands) mean_recourse = mean_recourse core_t_min = minimum(core_t) core_t_max = maximum(core_t)

# -----------------------------------------------------------------------------
# algorithm parameters
# -----------------------------------------------------------------------------
benders_param = BendersBnBParam(
    time_limit = time_limit,
    gap_tolerance = 1e-6,
    verbose = true,
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

# -----------------------------------------------------------------------------
# master model
# -----------------------------------------------------------------------------
master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)

# -----------------------------------------------------------------------------
# typical oracles: separable knapsack over scenarios (two copies for disjunction)
# -----------------------------------------------------------------------------
typical_oracles = [
    SeparableOracle(data, master, CFLKnapsackOracle(), data.n_scenarios; customize = customize_sub_model!, optimizer = optimizer),
    SeparableOracle(data, master, CFLKnapsackOracle(), data.n_scenarios; customize = customize_sub_model!, optimizer = optimizer),
]
disjunctive_oracle = DirectionalPolarDCGLPOracle(master, typical_oracles, oracle_param)

# -----------------------------------------------------------------------------
# lazy oracle + root preprocessing
# -----------------------------------------------------------------------------
lazy_oracle = SeparableOracle(data, master, CFLKnapsackOracle(), data.n_scenarios; customize = customize_sub_model!, optimizer = optimizer)
root_preprocessing = RootNodePreprocessing(
    lazy_oracle,
    BendersSeqInOut,
    BendersSeqInOutParam(
        time_limit = min(300.0, time_limit),
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
    @info "DirectionalPolarDCGLP SCFLP knapsack script (Gurobi, FLCAP) build completed without solve." instance = instance
else
    solve!(env)
    obj_value = try
        env.obj_value
    catch
        NaN
    end
    @info "DirectionalPolarDCGLP SCFLP knapsack script (Gurobi, FLCAP) finished" instance = instance termination_status = env.termination_status objective_value = obj_value
end
