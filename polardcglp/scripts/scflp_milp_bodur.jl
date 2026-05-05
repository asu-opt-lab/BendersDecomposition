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

# load settings
options, _ = parse_script_args(ARGS)

instance = get_string_option(options, "instance", "cap102-s250")
seed = get_int_option(options, "seed", 1)
output_dir = get_string_option(options, "output_dir", "output")
time_limit = get_float_option(options, "time_limit", 14400.0)

Random.seed!(seed)

const SCFLP_DIR = normpath(joinpath(@__DIR__, "..", "data", "SCFLP_bodur"))

@info "SCFLP MILP script (CPLEX, Bodur formulation)" instance = instance dataset_dir = SCFLP_DIR seed = seed time_limit = time_limit

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
data = read_scflp_bodur(instance; filepath = SCFLP_DIR)

# -----------------------------------------------------------------------------
# MIP model
# -----------------------------------------------------------------------------
mip_model = Model(mip_optimizer)
customize_mip_model!(mip_model, data)
set_time_limit_sec(mip_model, time_limit)
set_optimizer_attribute(mip_model, MOI.Silent(), false)
optimize!(mip_model)

@info termination_status(mip_model)
@info "Node count: $(node_count(mip_model))"
@info "Elapsed time: $(solve_time(mip_model))"
@info "Objective value: $(objective_value(mip_model))"
@info "Objective bound: $(objective_bound(mip_model))"
@info "Relative gap: $(relative_gap(mip_model))"
