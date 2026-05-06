using JuMP, DataFrames, Logging, CSV
using BendersX
using ArgParse
using Random
using Printf
using Statistics
using Gurobi

# Gurobi-based optimizers (replaces solver_defaults.jl which uses CPLEX)
optimizer = optimizer_with_attributes(
    Gurobi.Optimizer,
    "Threads" => 7,
    "FeasibilityTol" => 1e-9,
    "OptimalityTol" => 1e-9,
    "NumericFocus" => 1,
    MOI.Silent() => true,
)

mip_optimizer = optimizer_with_attributes(
    Gurobi.Optimizer,
    "Threads" => 7,
    "IntFeasTol" => 1e-9,
    "FeasibilityTol" => 1e-9,
    "MIPGap" => 1e-6,
    "OptimalityTol" => 1e-9,
    "NumericFocus" => 1,
    MOI.Silent() => true,
)

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--instance_no"
            help = "SNIP instance number"
            arg_type = Int
            default = 0
        "--snip_no"
            help = "SNIP interdiction setting"
            arg_type = Int
            default = 4
        "--budget"
            help = "Sensor installation budget"
            arg_type = Float64
            default = 30.0
        "--seed"
            help = "Random seed"
            arg_type = Int
            default = 1
        "--output_dir"
            help = "Output directory"
            arg_type = String
            default = "output"
    end
    return parse_args(s)
end

global_logger(ConsoleLogger(stderr, Logging.Debug))

time_limit = 14400.0

# load settings
args = parse_commandline()

Random.seed!(args["seed"])
instance_no = args["instance_no"]
snip_no = args["snip_no"]
budget = args["budget"]
output_dir = args["output_dir"]

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
data = read_snip_data(instance_no, snip_no, budget)

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
# classical oracle (separable across scenarios)
# -----------------------------------------------------------------------------
oracle = SeparableOracle(data, master, ClassicalOracle(), data.num_scenarios; customize = customize_sub_model!, optimizer = optimizer)

# -----------------------------------------------------------------------------
# root node preprocessing (BendersSeq)
# -----------------------------------------------------------------------------
root_seq_type = BendersSeq
root_param = BendersSeqParam(
    time_limit = 300.0,
    gap_tolerance = 1e-9,
    verbose = true
)
root_preprocessing = RootNodePreprocessing(oracle, root_seq_type, root_param)
# root_preprocessing = NoRootNodePreprocessing()

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
