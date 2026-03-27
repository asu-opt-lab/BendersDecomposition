using BendersDecomposition
using JuMP
using CPLEX
using Gurobi
using CuPDLPx
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
dcglp_solver_param = Dict("solver" => "CPLEX", "CPXPARAM_Threads" => 7, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1)
oracle_solver_param = Dict("solver" => "CPLEX", "CPXPARAM_Threads" => 7, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1)

# CuPDLPx
oracle_solver_param = Dict("solver" => "CuPDLPx", "verbose" => true, "iteration_limit" => 999999999)

# Gurobi PDHG
oracle_solver_param = Dict("solver" => "Gurobi", "method" => 6, "PDHGGPU" => 1, "LPWarmStart" => 0, "Crossover" => 1)


oracle_param = DisjunctiveOracleParam(
    norm = LpNorm(1.0), 
    split_index_selection_rule = RandomFractional(),
    disjunctive_cut_append_rule = AllDisjunctiveCuts(), 
    strengthened = true, 
    add_benders_cuts_to_master = true, 
    fraction_of_benders_cuts_to_master = 1.0, 
    reuse_dcglp = true,
    lift = false
)

dcglp_param = DcglpParam(
    time_limit = 200.0,
    gap_tolerance = 1e-3,
    halt_limit = 3,
    iter_limit = 15,
    verbose = true
)

user_cb_param = UserCallbackParam(frequency=500)
# -----------------------------------------------------------------------------
# autodecompose
# -----------------------------------------------------------------------------
data, master, disjunctive_oracle = auto_decompose(
    model,
    :disjunctive;
    master_solver_param = master_solver_param,
    oracle_param = oracle_param,
    typical_oracle_solver_param = oracle_solver_param,
    dcglp_solver_param = dcglp_solver_param,
    dcglp_param = dcglp_param
)

# -----------------------------------------------------------------------------
# lazy callback
# -----------------------------------------------------------------------------
_, _, typical_oracle = auto_decompose_unified(
    model;
    master_solver_param = master_solver_param,
    oracle_solver_param = oracle_solver_param
)

lazy_callback = LazyCallback(typical_oracle)

# -----------------------------------------------------------------------------
# user callback
# -----------------------------------------------------------------------------
user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)

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