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

instance_no = get_int_option(options, "instance_no", 0)
snip_no = get_int_option(options, "snip_no", 3)
budget = get_float_option(options, "budget", 30.0)
seed = get_int_option(options, "seed", 1)
output_dir = get_string_option(options, "output_dir", "output")
time_limit = get_float_option(options, "time_limit", 14400.0)
threads = get_int_option(options, "threads", 7)
build_only = get_bool_option(options, "build_only", false)

Random.seed!(seed)

sub_optimizer = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_Threads" => 1,
    "CPX_PARAM_EPRHS" => 1e-9,
    "CPX_PARAM_EPOPT" => 1e-9,
    "CPX_PARAM_NUMERICALEMPHASIS" => 1,
    MOI.Silent() => true,
)

@info "SNIP classical script (CPLEX)" instance_no = instance_no snip_no = snip_no budget = budget seed = seed time_limit = time_limit threads = threads build_only = build_only

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
data = read_snip_data(instance_no, snip_no, budget)

# -----------------------------------------------------------------------------
# master model
# -----------------------------------------------------------------------------
master = Master(data; customize = customize_master_model!, optimizer = master_optimizer)

# -----------------------------------------------------------------------------
# typical oracle: separable classical over scenarios
# -----------------------------------------------------------------------------
typical_oracle = SeparableOracle(
    data,
    master,
    ClassicalOracle(),
    data.num_scenarios;
    customize = customize_sub_model!,
    optimizer = sub_optimizer,
)

# -----------------------------------------------------------------------------
# root node preprocessing
# -----------------------------------------------------------------------------
root_preprocessing = RootNodePreprocessing(
    typical_oracle,
    BendersSeq,
    BendersSeqParam(
        time_limit = min(100.0, time_limit),
        gap_tolerance = 1e-6,
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
    @info "SNIP classical script (CPLEX) build completed without solve." instance_no = instance_no snip_no = snip_no budget = budget
else
    solve!(env)
    obj_value = try
        env.obj_value
    catch
        NaN
    end
    @info "SNIP classical script (CPLEX) finished" instance_no = instance_no snip_no = snip_no budget = budget termination_status = env.termination_status objective_value = obj_value
end
