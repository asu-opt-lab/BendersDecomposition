using JuMP, DataFrames, Logging, CSV
using BendersX
using Random
using Printf
using Statistics
using Gurobi

include(normpath(joinpath(@__DIR__, "..", "solver_defaults_gurobi.jl")))
include(normpath(joinpath(@__DIR__, "..", "script_utils.jl")))

global_logger(ConsoleLogger(stderr, Logging.Debug))

options, _ = parse_script_args(ARGS)

instance = get_string_option(options, "instance", "T100x100_5_1")
seed = get_int_option(options, "seed", 1)
output_dir = get_string_option(options, "output_dir", "output")
time_limit = get_float_option(options, "time_limit", 14400.0)
threads = get_int_option(options, "threads", 7)
build_only = get_bool_option(options, "build_only", false)

Random.seed!(seed)

@info "CFLP knapsack script (Gurobi)" instance = instance seed = seed time_limit = time_limit threads = threads build_only = build_only



# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
data = read_cfl_file(instance)

# -----------------------------------------------------------------------------
# master model
# -----------------------------------------------------------------------------
master = Master(data; customize = customize_master_model!, optimizer = master_optimizer)

# -----------------------------------------------------------------------------
# typical oracle: CFLKnapsack (deterministic, no scenario separation)
# -----------------------------------------------------------------------------
typical_oracle = CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = sub_optimizer)

# -----------------------------------------------------------------------------
# root node preprocessing
# -----------------------------------------------------------------------------
root_preprocessing = RootNodePreprocessing(
    typical_oracle,
    BendersSeqInOut,
    BendersSeqInOutParam(
        time_limit = min(300.0, time_limit),
        gap_tolerance = 1e-9,
        stabilizing_x = ones(data.n_facilities),
        α = 0.9,
        λ = 0.1,
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
    @info "CFLP knapsack script (Gurobi) build completed without solve." instance = instance
else
    solve!(env)
    obj_value = try
        env.obj_value
    catch
        NaN
    end
    @info "CFLP knapsack script (Gurobi) finished" instance = instance termination_status = env.termination_status objective_value = obj_value
end
