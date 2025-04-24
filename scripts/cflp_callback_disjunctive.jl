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
instance = "f100-c100-r3-1"
output_dir = "scripts"
# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
problem = read_GK_data(instance)
dim_x = problem.n_facilities
dim_t = 1
c_x = problem.fixed_costs
c_t = [1]
data = Data(dim_x, dim_t, problem, c_x, c_t)

# -----------------------------------------------------------------------------
# load parameters
# -----------------------------------------------------------------------------
# Algorithm parameters
benders_param = BendersBnBParam(
    time_limit = 200.0,
    gap_tolerance = 1e-6,
    verbose = true
)

dcglp_param = DcglpParam(
    time_limit = 1000.0,
    gap_tolerance = 1e-3,
    halt_limit = 3,
    iter_limit = 250,
    verbose = true
)

# Solver parameters
master_solver_param = Dict(
    "solver" => "CPLEX", 
    "CPX_PARAM_EPINT" => 1e-9, 
    "CPX_PARAM_EPRHS" => 1e-9
)

typical_oracle_solver_param = Dict(
    "solver" => "CPLEX", 
    "CPX_PARAM_EPRHS" => 1e-9, 
    "CPX_PARAM_NUMERICALEMPHASIS" => 1, 
    "CPX_PARAM_EPOPT" => 1e-9
)

dcglp_solver_param = Dict(
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
typical_oracles = [
    CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param),
    CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
]

for k = 1:2
    update_model!(typical_oracles[k], data)
end

# -----------------------------------------------------------------------------
# disjunctive oracle
# -----------------------------------------------------------------------------
disjunctive_oracle = DisjunctiveOracle(
    data, 
    typical_oracles; 
    solver_param = dcglp_solver_param, 
    param = dcglp_param
) 

oracle_param = DisjunctiveOracleParam(
    norm = LpNorm(1.0), 
    split_index_selection_rule = RandomFractional(),
    disjunctive_cut_append_rule = AllDisjunctiveCuts(), 
    strengthened = true, 
    add_benders_cuts_to_master = true, 
    fraction_of_benders_cuts_to_master = 0.5, 
    reuse_dcglp = true
)

set_parameter!(disjunctive_oracle, oracle_param)
update_model!(disjunctive_oracle, data)

# -----------------------------------------------------------------------------
# root node preprocessing
# -----------------------------------------------------------------------------
root_seq_type = BendersSeqInOut
root_param = BendersSeqInOutParam(
    time_limit = 100.0,
    gap_tolerance = 1e-6,
    stabilizing_x = ones(data.dim_x),
    α = 0.9,
    λ = 0.1,
    verbose = true
)

lazy_oracle = CFLKnapsackOracle(data; solver_param = typical_oracle_solver_param)
update_model!(lazy_oracle, data)

# Create root node preprocessing with oracle
root_preprocessing = RootNodePreprocessing(lazy_oracle, root_seq_type, root_param)

# -----------------------------------------------------------------------------
# lazy callback
# -----------------------------------------------------------------------------
lazy_callback = LazyCallback(lazy_oracle)

# -----------------------------------------------------------------------------
# user callback
# -----------------------------------------------------------------------------
user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=250))

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






