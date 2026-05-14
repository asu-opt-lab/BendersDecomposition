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

Random.seed!(seed)

@info "SNIP MILP script (CPLEX)" instance_no = instance_no snip_no = snip_no budget = budget seed = seed time_limit = time_limit threads = threads

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
data = read_snip_data(instance_no, snip_no, budget)

# -----------------------------------------------------------------------------
# MIP model
# -----------------------------------------------------------------------------
mip_model = Model(mip_optimizer)
customize_mip_model!(mip_model, data)
set_optimizer_attribute(mip_model, "CPXPARAM_Threads", threads)
set_time_limit_sec(mip_model, time_limit)
set_optimizer_attribute(mip_model, MOI.Silent(), false)
optimize!(mip_model)

@info termination_status(mip_model)
@info "Node count: $(node_count(mip_model))"
@info "Elapsed time: $(solve_time(mip_model))"
@info "Objective value: $(objective_value(mip_model))"
@info "Objective bound: $(objective_bound(mip_model))"
@info "Relative gap: $(relative_gap(mip_model))"
