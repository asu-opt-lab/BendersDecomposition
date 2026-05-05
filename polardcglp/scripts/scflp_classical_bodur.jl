using JuMP, DataFrames, Logging, CSV
using BendersX
using Random
using Printf
using Statistics
using CPLEX

include(normpath(joinpath(@__DIR__, "solver_defaults.jl")))
include(normpath(joinpath(@__DIR__, "script_utils.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "newscflp_data.jl")))
include(normpath(joinpath(@__DIR__, "..", "src", "newscflp_model.jl")))

global_logger(ConsoleLogger(stderr, Logging.Debug))

options, _ = parse_script_args(ARGS)

instance = get_string_option(options, "instance", "cap102-s250")
seed = get_int_option(options, "seed", 1)
output_dir = get_string_option(options, "output_dir", "output")
time_limit = get_float_option(options, "time_limit", 14400.0)
build_only = get_bool_option(options, "build_only", false)

Random.seed!(seed)

const SCFLP_DIR = normpath(joinpath(@__DIR__, "..", "data", "SCFLP_bodur"))

@info "SCFLP classical script (CPLEX, Bodur formulation)" instance = instance dataset_dir = SCFLP_DIR seed = seed time_limit = time_limit build_only = build_only

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
data = read_scflp_bodur(instance; filepath = SCFLP_DIR)

# -----------------------------------------------------------------------------
# master model
# -----------------------------------------------------------------------------
master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)

# -----------------------------------------------------------------------------
# typical oracle: separable classical over scenarios
# -----------------------------------------------------------------------------
typical_oracle = SeparableOracle(
    data,
    master,
    ClassicalOracle(),
    data.n_scenarios;
    customize = customize_sub_model!,
    optimizer = optimizer,
)

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
    @info "SCFLP classical script (CPLEX, Bodur formulation) build completed without solve." instance = instance
else
    solve!(env)
    obj_value = try
        env.obj_value
    catch
        NaN
    end
    @info "SCFLP classical script (CPLEX, Bodur formulation) finished" instance = instance termination_status = env.termination_status objective_value = obj_value
end
