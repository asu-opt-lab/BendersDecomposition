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
        "--instance_no"
            help = "SNIP instance number"
            arg_type = Int
            default = 0
        "--snip_no"
            help = "SNIP interdiction setting"
            arg_type = Int
            default = 3
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

threads = 7
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
# MIP model
# -----------------------------------------------------------------------------
mip_model = Model(mip_optimizer)
customize_mip_model!(mip_model, data)
set_optimizer_attribute(mip_model, "Threads", threads)
set_time_limit_sec(mip_model, time_limit)
set_optimizer_attribute(mip_model, MOI.Silent(), false)
optimize!(mip_model)

@info termination_status(mip_model)
@info "Node count: $(node_count(mip_model))"
@info "Elapsed time: $(solve_time(mip_model))"
@info "Objective value: $(objective_value(mip_model))"
@info "Objective bound: $(objective_bound(mip_model))"
@info "Relative gap: $(relative_gap(mip_model))"
