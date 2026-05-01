using JuMP, DataFrames, Logging, CSV
using BendersX
using ArgParse
using Random
using Printf
using Statistics
using Gurobi

# Gurobi-based mip_optimizer (replaces solver_defaults.jl which uses CPLEX)
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
    optimizer = optimizer_with_attributes(Gurobi.Optimizer,
        "Threads" => 1, "IntFeasTol" => 1e-9,
        "FeasibilityTol" => 1e-9, "MIPGap" => 1e-6,
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
    time_limit = 14000.0,
    gap_tolerance = 1e-6,
    verbose = true
)

# -----------------------------------------------------------------------------
# master model
# -----------------------------------------------------------------------------
master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)
set_optimizer_attribute(master.model, "BranchDir", 1)

# -----------------------------------------------------------------------------
# typical oracle
# -----------------------------------------------------------------------------
typical_oracle = UFLKnapsackOracle(data)
set_parameter!(typical_oracle, "add_only_violated_cuts", true)

# -----------------------------------------------------------------------------
# root node preprocessing
# -----------------------------------------------------------------------------
root_seq_type = BendersSeq
root_param = BendersSeqParam(
    time_limit = 100.0,
    gap_tolerance = 1e-9,
    verbose = true
)

# Create root node preprocessing with oracle
root_preprocessing = RootNodePreprocessing(typical_oracle, root_seq_type, root_param)

# -----------------------------------------------------------------------------
# lazy callback
# -----------------------------------------------------------------------------
lazy_callback = LazyCallback(typical_oracle)

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
