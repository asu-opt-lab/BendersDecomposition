using JuMP, DataFrames, Logging, CSV
using BendersX
using ArgParse
using Printf
using Statistics
using CPLEX
include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--instance"
            help = "Instance name"
            arg_type = String
            default = "T100x100_5_1"
        "--time_limit"
            help = "Overall Benders BnB time limit in seconds"
            arg_type = Float64
            default = 3600.0
        "--root_time_limit"
            help = "Root-node classical Benders preprocessing time limit in seconds"
            arg_type = Float64
            default = 3600.0
        "--output_dir"
            help = "Output directory"
            arg_type = String
            default = "output"
    end
    return parse_args(s)
end

global_logger(ConsoleLogger(stderr, Logging.Debug))

# load settings
args = parse_commandline()

instance = args["instance"]
time_limit = args["time_limit"]
root_time_limit = args["root_time_limit"]
output_dir = args["output_dir"]

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
data = read_cfl_file(instance)

# -----------------------------------------------------------------------------
# load parameters
# -----------------------------------------------------------------------------
benders_param = BendersBnBParam(
    time_limit = time_limit,
    gap_tolerance = 1e-6,
    verbose = true
)

# -----------------------------------------------------------------------------
# master model
# -----------------------------------------------------------------------------
master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)

# -----------------------------------------------------------------------------
# classical oracle
# -----------------------------------------------------------------------------
oracle = ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer)

# -----------------------------------------------------------------------------
# root node preprocessing (classical sequential Benders)
# -----------------------------------------------------------------------------
root_seq_type = BendersSeq
root_param = BendersSeqParam(
    time_limit = root_time_limit,
    gap_tolerance = 1e-9,
    verbose = true
)
root_preprocessing = RootNodePreprocessing(oracle, root_seq_type, root_param)

# -----------------------------------------------------------------------------
# lazy callback
# -----------------------------------------------------------------------------
lazy_callback = LazyCallback(oracle)

# -----------------------------------------------------------------------------
# user callback
# -----------------------------------------------------------------------------
user_callback = NoUserCallback()

# -----------------------------------------------------------------------------
# BendersBnB
# -----------------------------------------------------------------------------
env = BendersBnB(
    master,
    root_preprocessing,
    lazy_callback,
    user_callback;
    param = benders_param
)

# -----------------------------------------------------------------------------
# solve
# -----------------------------------------------------------------------------
solution_log = solve!(env)

@info env.termination_status
@info "Objective value: $(env.obj_value)"
