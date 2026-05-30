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
# master model
# -----------------------------------------------------------------------------
master = Master(data; model = update_master_model!, optimizer = mip_optimizer)

# -----------------------------------------------------------------------------
# typical oracle
# -----------------------------------------------------------------------------
typical_oracle = CFLKnapsackOracle(data, master; model = update_sub_model!, optimizer = optimizer)

# -----------------------------------------------------------------------------
# BendersSeqInOut
# -----------------------------------------------------------------------------
start_time = time()
undo = relax_integrality(master.model)
benders_inout_param = BendersSeqInOutParam(;
            time_limit = 300.0,
            gap_tolerance = 1e-6,
            verbose = true,
            stabilizing_x = ones(data.n_facilities),
            α = 0.9,
            λ = 0.1
            )
set_optimizer_attribute(typical_oracle.model, "CPX_PARAM_LPMETHOD", 2)
env = BendersSeqInOut(master, typical_oracle; param = benders_inout_param)
solution_log = solve!(env)
set_optimizer_attribute(typical_oracle.model, "CPX_PARAM_LPMETHOD", 0)
println("BendersSeqInOut done")
undo()
spend_time_prev = time() - start_time
println("Spend time: $spend_time_prev seconds")

# -----------------------------------------------------------------------------
# BendersBnB
# -----------------------------------------------------------------------------
root_preprocessing = NoRootNodePreprocessing()
lazy_callback = LazyCallback(typical_oracle)
user_callback = NoUserCallback()

benders_param = BendersBnBParam(
    time_limit = 14400.0 - spend_time_prev,
    gap_tolerance = 1e-6,
    verbose = true
)
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
println("BendersBnB done")
@info "total time: $(time() - start_time) seconds"
