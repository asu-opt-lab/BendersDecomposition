using JuMP, DataFrames, Logging, CSV
using BendersDecomposition
using Printf  
using Statistics  
import BendersDecomposition: generate_cuts
include("$(dirname(@__DIR__))/example/cflp/data_reader.jl")
include("$(dirname(@__DIR__))/example/cflp/oracle.jl")
include("$(dirname(@__DIR__))/example/cflp/model.jl")


# load settings
# args = parse_commandline()

# instance = args["instance"]
# output_dir = args["output_dir"]
instance = "T100x100_5_1"
output_dir = "scripts"
# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
problem = read_cfl_file(instance)
dim_x = problem.n_facilities
dim_t = 1
c_x = problem.fixed_costs
c_t = [1]
data = Data(dim_x, dim_t, problem, c_x, c_t)

# -----------------------------------------------------------------------------
# load parameters
# -----------------------------------------------------------------------------
# Algorithm parameters

# Solver parameters
master_solver_param = Dict(
    "solver" => "CPLEX", 
    "CPX_PARAM_EPINT" => 1e-9, 
    "CPX_PARAM_EPRHS" => 1e-9,
    "CPX_PARAM_EPGAP" => 1e-6
)

typical_oracle_solver_param = Dict(
    "solver" => "CPLEX", 
    "CPX_PARAM_EPRHS" => 1e-9, 
    "CPX_PARAM_NUMERICALEMPHASIS" => 1, 
    "CPX_PARAM_EPOPT" => 1e-9
)

# -----------------------------------------------------------------------------
# master model
# -----------------------------------------------------------------------------
master = Master(data; solver_param = master_solver_param)
update_model!(master, data)

# -----------------------------------------------------------------------------
# typical oracles
# -----------------------------------------------------------------------------
# Create two oracles for kappa & nu
typical_oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
update_model!(typical_oracle, data)

# -----------------------------------------------------------------------------
# lazy callback
# -----------------------------------------------------------------------------
lazy_callback = LazyCallback(typical_oracle)

# -----------------------------------------------------------------------------
# user callback
# -----------------------------------------------------------------------------
user_callback = NoUserCallback()

# -----------------------------------------------------------------------------
# BendersSeqInOut
# -----------------------------------------------------------------------------
start_time = time()
undo = relax_integrality(master.model)
benders_inout_param = BendersSeqInOutParam(;
            time_limit = 300.0,
            gap_tolerance = 1e-6,
            verbose = true,
            stabilizing_x = ones(data.dim_x),
            α = 0.9,
            λ = 0.1
            )
set_optimizer_attribute(typical_oracle.model, "CPX_PARAM_LPMETHOD", 2) 
env = BendersSeqInOut(data, master, typical_oracle; param = benders_inout_param)
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
benders_param = BendersBnBParam(
    time_limit = 14400.0 - spend_time_prev,
    gap_tolerance = 1e-6,
    verbose = true
)
env = BendersBnB(
    data, 
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





