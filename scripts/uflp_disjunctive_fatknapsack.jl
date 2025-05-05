using JuMP, DataFrames, Logging, CSV
using BendersDecomposition
using Random
using Printf  
using Statistics  
import BendersDecomposition: generate_cuts
include("$(dirname(@__DIR__))/example/uflp/data_reader.jl")
include("$(dirname(@__DIR__))/example/uflp/oracle.jl")
include("$(dirname(@__DIR__))/example/uflp/model.jl")
Random.seed!(1234)

# load settings
args = parse_commandline()

instance = args["instance"]
output_dir = args["output_dir"]

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
problem = read_Simple_data(instance)
dim_x = problem.n_facilities
dim_t = problem.n_customers
c_x = problem.fixed_costs
c_t = ones(dim_t)
data = Data(dim_x, dim_t, problem, c_x, c_t)

# -----------------------------------------------------------------------------
# load parameters
# -----------------------------------------------------------------------------
# Algorithm parameters
benders_param = BendersBnBParam(
    time_limit = 3600.0,
    gap_tolerance = 1e-6,
    #disjunctive_root_process = true,
    verbose = true
)

dcglp_param = DcglpParam(
    time_limit = 1000.0,
    gap_tolerance = 1e-3,
    halt_limit = 3,
    iter_limit = 3,
    verbose = true
)

# Solver parameters
master_solver_param = Dict(
    "solver" => "CPLEX", 
    "CPX_PARAM_EPINT" => 1e-9, 
    "CPX_PARAM_EPRHS" => 1e-9,
    "CPX_PARAM_EPGAP" => 1e-6
)

dcglp_solver_param = Dict(
    "solver" => "CPLEX", 
    "CPX_PARAM_EPRHS" => 1e-9, 
    "CPX_PARAM_NUMERICALEMPHASIS" => 1, 
    "CPX_PARAM_EPOPT" => 1e-9
) #LPMETHOD default 0, network simplex 3

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
    UFLKnapsackOracle(data),
    UFLKnapsackOracle(data)
]

for k=1:2
    set_parameter!(typical_oracles[k], "add_only_violated_cuts", true)
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
    norm = LpNorm(Inf), 
    split_index_selection_rule = RandomFractional(),
    disjunctive_cut_append_rule = AllDisjunctiveCuts(), 
    strengthened = true, 
    add_benders_cuts_to_master = true,  
    fraction_of_benders_cuts_to_master = 0.05, 
    reuse_dcglp = false, 
    lift = true
) 

set_parameter!(disjunctive_oracle, oracle_param)
update_model!(disjunctive_oracle, data)

# -----------------------------------------------------------------------------
# root node preprocessing
# -----------------------------------------------------------------------------
root_seq_type = BendersSeq
root_param = BendersSeqParam(
    time_limit = 100.0,
    gap_tolerance = 1e-6,
    verbose = true
)

lazy_oracle = UFLKnapsackOracle(data)
set_parameter!(lazy_oracle, "add_only_violated_cuts", true)

# Create root node preprocessing with oracle
root_preprocessing = RootNodePreprocessing(lazy_oracle, root_seq_type, root_param)

# -----------------------------------------------------------------------------
# lazy callback
# -----------------------------------------------------------------------------
lazy_callback = LazyCallback(lazy_oracle)

# -----------------------------------------------------------------------------
# user callback
# -----------------------------------------------------------------------------
user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=500))

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






