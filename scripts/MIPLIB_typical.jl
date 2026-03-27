using BendersDecomposition
using JuMP
using CPLEX
using CuPDLPx
using Gurobi
using ArgParse
using Random
Random.seed!(1234)

# -----------------------------------------------------------------------------
# command line parsing: ONLY --instance
# -----------------------------------------------------------------------------
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--instance"
            help = "Instance name (e.g., r01.1.dow)"
            arg_type = String
            required = true
    end

    return parse_args(s)
end

# load settings
args = parse_commandline()
instance = args["instance"]

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
path = joinpath("scripts", "MIPLIB_data", instance, instance * ".mps.gz")
model = read_from_file(path)

# -----------------------------------------------------------------------------
# load parameters
# -----------------------------------------------------------------------------
# CPLEX
master_solver_param = Dict("solver" => "CPLEX", "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6)
oracle_solver_param = Dict("solver" => "CPLEX", "CPXPARAM_Threads" => 7, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1)

# PDHG
oracle_solver_param = Dict("solver" => "Gurobi", "method" => 6, "PDHGGPU" => 1, "LPWarmStart" => 0, "Threads" => 7, "Crossover" => 1)


# -----------------------------------------------------------------------------
# autodecompose
# -----------------------------------------------------------------------------
data, master, typical_oracle = auto_decompose_unified(model; master_solver_param = master_solver_param, oracle_solver_param = oracle_solver_param)

# -----------------------------------------------------------------------------
# lazy callback
# -----------------------------------------------------------------------------
lazy_callback = LazyCallback(typical_oracle)

# -----------------------------------------------------------------------------
# user callback
# -----------------------------------------------------------------------------
user_callback = NoUserCallback()

# -----------------------------------------------------------------------------
# BendersSeq
# -----------------------------------------------------------------------------
start_time = time()
undo = relax_integrality(master.model)

root_param = BendersSeqParam(;
    time_limit = 300.0,
    gap_tolerance = 1e-6,
    verbose = true
)

env = BendersSeq(data, master, typical_oracle; param = root_param)
solution_log = solve!(env)

undo()
spend_time_prev = time() - start_time
println("Spend time: $spend_time_prev seconds")

# -----------------------------------------------------------------------------
# BendersBnB
# -----------------------------------------------------------------------------
root_preprocessing = NoRootNodePreprocessing()
benders_param = BendersBnBParam(
    time_limit = 3600.0 - spend_time_prev,
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
println("termination status: $(termination_status(env.master.model))")
println("total time: $(time() - start_time) seconds")