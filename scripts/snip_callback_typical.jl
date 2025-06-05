using JuMP, DataFrames, Logging, CSV
using BendersDecomposition
using Printf  
using Statistics  
import BendersDecomposition: generate_cuts
include("$(dirname(@__DIR__))/example/snip/data_reader.jl")
include("$(dirname(@__DIR__))/example/snip/model.jl")


# load settings
args = parse_commandline()

instance = args["snip_instance"]
snipno = args["snip_no"]
budget = args["snip_budget"]
println("instance: $instance, snipno: $snipno, budget: $budget")

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
problem = read_snip_data(instance, snipno, budget)
dim_x = length(problem.D)
dim_t = problem.num_scenarios
c_x = zeros(dim_x)
c_t = map(p -> p[3], problem.scenarios)
data = Data(dim_x, dim_t, problem, c_x, c_t)

# -----------------------------------------------------------------------------
# load parameters
# -----------------------------------------------------------------------------
# Algorithm parameters
benders_param = BendersBnBParam(
    time_limit = 3600.0,
    gap_tolerance = 1e-6,
    verbose = true
)

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
    "CPX_PARAM_EPOPT" => 1e-9,
    "CPX_PARAM_THREADS" => 7
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
typical_oracle = SeparableOracle(data, ClassicalOracle(), data.problem.num_scenarios; solver_param = typical_oracle_solver_param)
for j=1:typical_oracle.N
    update_model!(typical_oracle.oracles[j], data, j)
end

# -----------------------------------------------------------------------------
# root node preprocessing
# -----------------------------------------------------------------------------
# root_seq_type = BendersSeqInOut
# root_param = BendersSeqInOutParam(
#     time_limit = 300.0,
#     gap_tolerance = 1e-6,
#     stabilizing_x = ones(data.dim_x),
#     α = 0.9,
#     λ = 0.1,
#     verbose = true
# )

root_seq_type = BendersSeq
root_param = BendersSeqParam(
    time_limit = 300.0,
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






