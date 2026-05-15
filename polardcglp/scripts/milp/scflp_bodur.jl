using JuMP, DataFrames, Logging, CSV
using BendersX
using Random
using Printf
using Statistics
using CPLEX

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))
include(normpath(joinpath(@__DIR__, "..", "script_utils.jl")))

global_logger(ConsoleLogger(stderr, Logging.Debug))

# load settings
options, _ = parse_script_args(ARGS)

instance = get_string_option(options, "instance", "cap102-s250")
seed = get_int_option(options, "seed", 1)
output_dir = get_string_option(options, "output_dir", "output")
time_limit = get_float_option(options, "time_limit", 14400.0)
threads = get_int_option(options, "threads", 7)
build_only = get_bool_option(options, "build_only", false)

Random.seed!(seed)

const SCFLP_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "SCFLP_bodur"))

@info "SCFLP MILP script (CPLEX, Bodur formulation)" instance = instance dataset_dir = SCFLP_DIR seed = seed time_limit = time_limit threads = threads build_only = build_only

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
data = read_scflp_bodur(instance; filepath = SCFLP_DIR)

# -----------------------------------------------------------------------------
# MIP model
# -----------------------------------------------------------------------------
mip_optimizer_local = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_Threads" => threads,
    "CPX_PARAM_EPINT" => 1e-9,
    "CPX_PARAM_EPRHS" => 1e-9,
    "CPX_PARAM_EPGAP" => 1e-6,
    MOI.Silent() => true,
)

mip_model = Model(mip_optimizer_local)
customize_mip_model!(mip_model, data)
set_time_limit_sec(mip_model, time_limit)
set_optimizer_attribute(mip_model, MOI.Silent(), false)

if build_only
    @info "SCFLP MILP script build completed without solve." instance = instance
else
    optimize!(mip_model)

    @info termination_status(mip_model)
    @info "Node count: $(node_count(mip_model))"
    @info "Elapsed time: $(solve_time(mip_model))"
    @info "Objective value: $(objective_value(mip_model))"
    @info "Objective bound: $(objective_bound(mip_model))"
    @info "Relative gap: $(relative_gap(mip_model))"
end
