using JuMP, DataFrames, Logging, CSV
using BendersX
using Printf
using Statistics
using CPLEX
include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))

# load settings
# args = parse_commandline()

# instance = args["instance"]
# output_dir = args["output_dir"]
instance = "T100x100_5_1"
output_dir = "scripts"

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
data = read_cfl_file(instance)

# -----------------------------------------------------------------------------
# load parameters
# -----------------------------------------------------------------------------
# Algorithm parameters
benders_param = BendersBnBParam(
    time_limit = 200.0,
    gap_tolerance = 1e-6,
    verbose = true
)

dcglp_optimizer = optimizer_with_attributes(CPLEX.Optimizer,
    "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1,
    "CPX_PARAM_EPOPT" => 1e-9, MOI.Silent() => true)
dcglp_param = DcglpParam(dcglp_optimizer;
    time_limit = 1000.0,
    gap_tolerance = 1e-3,
    halt_limit = 3,
    iter_limit = 250,
    verbose = true
)

oracle_param = SplitOracleParam(LpDistanceNormalization(LpNorm(1.0)); dcglp_param = dcglp_param,
    split_index_selection_rule = RandomFractional(),
    disjunctive_cut_append_rule = AllDisjunctiveCuts(),
    strengthened = true,
    add_benders_cuts_to_master = true,
    fraction_of_benders_cuts_to_master = 0.5,
    reuse_dcglp = true
)

# -----------------------------------------------------------------------------
# master model
# -----------------------------------------------------------------------------
master = Master(data; model = update_master_model!, optimizer = mip_optimizer)

# -----------------------------------------------------------------------------
# typical oracles
# -----------------------------------------------------------------------------
# Create two oracles for kappa & nu
typical_oracles = [
    CFLKnapsackOracle(data, master; model = update_sub_model!, optimizer = optimizer),
    CFLKnapsackOracle(data, master; model = update_sub_model!, optimizer = optimizer)
]

# -----------------------------------------------------------------------------
# disjunctive oracle
# -----------------------------------------------------------------------------
disjunctive_oracle = SplitOracle(master, Tuple(typical_oracles), oracle_param)

# -----------------------------------------------------------------------------
# root node preprocessing
# -----------------------------------------------------------------------------
lazy_oracle = CFLKnapsackOracle(data, master; model = update_sub_model!, optimizer = optimizer)

root_seq_type = BendersSeqInOut
root_param = BendersSeqInOutParam(
    time_limit = 100.0,
    gap_tolerance = 1e-6,
    stabilizing_x = ones(data.n_facilities),
    α = 0.9,
    λ = 0.1,
    verbose = true
)

# Create root node preprocessing with oracle
root_preprocessing = RootNodePreprocessing(lazy_oracle, root_seq_type, root_param)

# -----------------------------------------------------------------------------
# lazy callback
# -----------------------------------------------------------------------------
lazy_callback = LazyCallback(lazy_oracle)

# -----------------------------------------------------------------------------
# user callback
# -----------------------------------------------------------------------------
user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=250))

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
