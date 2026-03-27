using JuMP, DataFrames, Logging, CSV
using BendersX
using ArgParse
using Random
using Printf
using Statistics
using CPLEX

global_logger(ConsoleLogger(stderr, Logging.Debug))

frequency = 500
reuse_dcglp = false
add_benders_cuts_to_master = 2
fraction_of_benders_cuts_to_master = 0.05
threads = 7
split_index_selection_rule = LargestFractional()
time_limit = 14400.0
p = Inf
adjust_t_to_fx = false
lift = false
strengthened = true

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--instance"
            help = "Instance name"
            arg_type = String
            default = "ga250a-1"
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

# load settings
args = parse_commandline()

Random.seed!(args["seed"])
instance = args["instance"]
output_dir = args["output_dir"]

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
data = read_Simple_data(instance)

# -----------------------------------------------------------------------------
# Customize master model function for knapsack oracle (needs t[1:n_customers])
# -----------------------------------------------------------------------------
function customize_master_model!(model::Model, data::UFLPData)
    optimizer = optimizer_with_attributes(CPLEX.Optimizer,
        "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9,
        "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6,
        MOI.Silent() => true)
    set_optimizer(model, optimizer)
    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t[1:data.n_customers] >= -1e6)
    @constraint(model, sum(x) >= 2)
    @objective(model, Min, data.fixed_costs'* x + sum(t))
    return (x = x, ), t
end

# -----------------------------------------------------------------------------
# load parameters
# -----------------------------------------------------------------------------
# Algorithm parameters
benders_param = BendersBnBParam(
    time_limit = time_limit,
    gap_tolerance = 1e-6,
    verbose = true
)

dcglp_optimizer = optimizer_with_attributes(CPLEX.Optimizer,
    "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1,
    "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_THREADS" => threads,
    MOI.Silent() => true)
dcglp_param = DcglpParam(dcglp_optimizer;
    time_limit = 100.0,
    gap_tolerance = 1e-3,
    halt_limit = 3,
    iter_limit = 250,
    verbose = true
)

oracle_param = SplitOracleParam(dcglp_param;
    norm = LpNorm(p),
    split_index_selection_rule = split_index_selection_rule,
    disjunctive_cut_append_rule = AllDisjunctiveCuts(),
    strengthened = true,
    add_benders_cuts_to_master = add_benders_cuts_to_master,
    fraction_of_benders_cuts_to_master = fraction_of_benders_cuts_to_master,
    reuse_dcglp = reuse_dcglp,
    lift = lift,
    adjust_t_to_fx = adjust_t_to_fx
)

# -----------------------------------------------------------------------------
# master model
# -----------------------------------------------------------------------------
master = Master(data; customize = customize_master_model!)
set_optimizer_attribute(master.model, "CPX_PARAM_BRDIR", 1)

# -----------------------------------------------------------------------------
# typical oracles
# -----------------------------------------------------------------------------
# Create two oracles for kappa & nu
typical_oracles = [
    UFLKnapsackOracle(data),
    UFLKnapsackOracle(data)
]

for k=1:2
    set_parameter!(typical_oracles[k], "add_only_violated_cuts", true)
end

# -----------------------------------------------------------------------------
# disjunctive oracle
# -----------------------------------------------------------------------------
disjunctive_oracle = SplitOracle(master, typical_oracles, oracle_param)

# -----------------------------------------------------------------------------
# root node preprocessing
# -----------------------------------------------------------------------------
lazy_oracle = UFLKnapsackOracle(data)
set_parameter!(lazy_oracle, "add_only_violated_cuts", true)

root_seq_type = BendersSeq
root_param = BendersSeqParam(
    time_limit = 100.0,
    gap_tolerance = 1e-9,
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
user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=frequency))

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
