using JuMP, DataFrames, Logging, CSV
using BendersX
using ArgParse
using Random
using Printf
using Statistics
using CPLEX
include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))

global_logger(ConsoleLogger(stderr, Logging.Debug))

frequency = 500
reuse_dcglp = false
add_benders_cuts_to_master = 2
fraction_of_benders_cuts_to_master = 0.05
threads = 7
split_index_selection_rule = LargestFractional()
time_limit = 14400.0
root_time_limit = 100.0
dcglp_time_limit = 100.0
dcglp_iter_limit = 250
dcglp_halt_limit = 3
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
        "--time_limit"
            help = "Benders branch-and-bound time limit in seconds"
            arg_type = Float64
            default = time_limit
        "--root_time_limit"
            help = "Root preprocessing time limit in seconds"
            arg_type = Float64
            default = root_time_limit
        "--dcglp_time_limit"
            help = "DCGLP time limit in seconds"
            arg_type = Float64
            default = dcglp_time_limit
        "--dcglp_iter_limit"
            help = "DCGLP iteration limit"
            arg_type = Int
            default = dcglp_iter_limit
        "--dcglp_halt_limit"
            help = "DCGLP halt limit"
            arg_type = Int
            default = dcglp_halt_limit
        "--frequency"
            help = "User callback frequency"
            arg_type = Int
            default = frequency
        "--threads"
            help = "Number of CPLEX threads"
            arg_type = Int
            default = threads
        "--reuse_dcglp"
            help = "Reuse DCGLP across callback calls"
            arg_type = Bool
            default = reuse_dcglp
        "--strengthened"
            help = "Use strengthened split cuts"
            arg_type = Bool
            default = strengthened
        "--lift"
            help = "Use lifted split cuts"
            arg_type = Bool
            default = lift
        "--adjust_t_to_fx"
            help = "Adjust t to f(x) before DCGLP separation"
            arg_type = Bool
            default = adjust_t_to_fx
        "--build_only"
            help = "Build the Benders environment without solving"
            arg_type = Bool
            default = false
    end
    return parse_args(s)
end

# load settings
args = parse_commandline()

Random.seed!(args["seed"])
instance = args["instance"]
output_dir = args["output_dir"]
time_limit = args["time_limit"]
root_time_limit = args["root_time_limit"]
dcglp_time_limit = args["dcglp_time_limit"]
dcglp_iter_limit = args["dcglp_iter_limit"]
dcglp_halt_limit = args["dcglp_halt_limit"]
frequency = args["frequency"]
threads = args["threads"]
reuse_dcglp = args["reuse_dcglp"]
strengthened = args["strengthened"]
lift = args["lift"]
adjust_t_to_fx = args["adjust_t_to_fx"]
build_only = args["build_only"]

@info "UFLP Benders BnB disjunctive knapsack script" instance = instance seed = args["seed"] time_limit = time_limit dcglp_time_limit = dcglp_time_limit frequency = frequency threads = threads reuse_dcglp = reuse_dcglp strengthened = strengthened lift = lift adjust_t_to_fx = adjust_t_to_fx build_only = build_only

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
data = read_Simple_data(instance)

# -----------------------------------------------------------------------------
# Customize master model function for knapsack oracle (needs t[1:n_customers])
# -----------------------------------------------------------------------------
function update_master_model!(model::Model, data::UFLPData)
    optimizer = optimizer_with_attributes(CPLEX.Optimizer,
        "CPXPARAM_Threads" => threads, "CPX_PARAM_EPINT" => 1e-9,
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
    time_limit = dcglp_time_limit,
    gap_tolerance = 1e-3,
    halt_limit = dcglp_halt_limit,
    iter_limit = dcglp_iter_limit,
    verbose = true
)

oracle_param = SplitOracleParam(LpDistanceNormalization(p); dcglp_param = dcglp_param,
    split_index_selection_rule = split_index_selection_rule,
    disjunctive_cut_append_rule = AllDisjunctiveCuts(),
    strengthened = strengthened,
    add_benders_cuts_to_master = add_benders_cuts_to_master,
    fraction_of_benders_cuts_to_master = fraction_of_benders_cuts_to_master,
    reuse_dcglp = reuse_dcglp,
    lift = lift,
    adjust_t_to_fx = adjust_t_to_fx
)

# -----------------------------------------------------------------------------
# master model
# -----------------------------------------------------------------------------
master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
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
disjunctive_oracle = SplitOracle(master, Tuple(typical_oracles); param = oracle_param)

# -----------------------------------------------------------------------------
# root node preprocessing
# -----------------------------------------------------------------------------
lazy_oracle = UFLKnapsackOracle(data)
set_parameter!(lazy_oracle, "add_only_violated_cuts", true)

root_seq_type = BendersSeq
root_param = BendersSeqParam(
    time_limit = min(root_time_limit, time_limit),
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
if build_only
    @info "UFLP Benders BnB disjunctive knapsack script build completed without solve." instance = instance
else
    solution_log = solve!(env)
    obj_value = try
        env.obj_value
    catch
        NaN
    end
    @info "UFLP Benders BnB disjunctive knapsack script finished" instance = instance termination_status = env.termination_status objective_value = obj_value
end
