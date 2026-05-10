using JuMP, DataFrames, Logging, CSV
using BendersX
using Random
using Printf
using Statistics
using CPLEX

include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))
include(normpath(joinpath(@__DIR__, "..", "script_utils.jl")))

global_logger(ConsoleLogger(stderr, Logging.Debug))

options, _ = parse_script_args(ARGS)

instance = get_string_option(options, "instance", "ga250a-2")
seed = get_int_option(options, "seed", 1)
output_dir = get_string_option(options, "output_dir", "output")
time_limit = get_float_option(options, "time_limit", 14400.0)
threads = get_int_option(options, "threads", 7)
build_only = get_bool_option(options, "build_only", false)

Random.seed!(seed)

@info "UFLP knapsack script (CPLEX)" instance = instance seed = seed time_limit = time_limit threads = threads build_only = build_only


# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
data = read_Simple_data(instance)

# -----------------------------------------------------------------------------
# master model (knapsack-style with sum(x) >= 2)
# -----------------------------------------------------------------------------
function customize_master_knapsack!(model::Model, data::UFLPData)
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer,
        "CPXPARAM_Threads" => threads,
        "CPX_PARAM_EPINT" => 1e-9,
        "CPX_PARAM_EPRHS" => 1e-9,
        "CPX_PARAM_EPGAP" => 1e-6,
        MOI.Silent() => true,
    )
    set_optimizer(model, optimizer)
    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t[1:data.n_customers] >= -1e6)
    @constraint(model, sum(x) >= 2)
    @objective(model, Min, data.fixed_costs' * x + sum(t))
    return (x = x,), t
end

master = Master(data; customize = customize_master_knapsack!, optimizer = master_optimizer)
set_optimizer_attribute(master.model, "CPX_PARAM_BRDIR", 1)

# -----------------------------------------------------------------------------
# typical oracle: UFLKnapsack (deterministic, no scenario separation)
# -----------------------------------------------------------------------------
typical_oracle = UFLKnapsackOracle(data)
set_parameter!(typical_oracle, "add_only_violated_cuts", true)

# -----------------------------------------------------------------------------
# root node preprocessing
# -----------------------------------------------------------------------------
root_preprocessing = RootNodePreprocessing(
    typical_oracle,
    BendersSeq,
    BendersSeqParam(
        time_limit = min(100.0, time_limit),
        gap_tolerance = 1e-9,
        verbose = true,
    ),
)

# -----------------------------------------------------------------------------
# callbacks
# -----------------------------------------------------------------------------
lazy_callback = LazyCallback(typical_oracle)
user_callback = NoUserCallback()

# -----------------------------------------------------------------------------
# BendersBnB
# -----------------------------------------------------------------------------
benders_param = BendersBnBParam(
    time_limit = time_limit,
    gap_tolerance = 1e-6,
    verbose = true,
)

env = BendersBnB(
    master,
    root_preprocessing,
    lazy_callback,
    user_callback;
    param = benders_param,
)

if build_only
    @info "UFLP knapsack script (CPLEX) build completed without solve." instance = instance
else
    solve!(env)
    obj_value = try
        env.obj_value
    catch
        NaN
    end
    @info "UFLP knapsack script (CPLEX) finished" instance = instance termination_status = env.termination_status objective_value = obj_value
end
