using JuMP, DataFrames, Logging, CSV
using BendersDecomposition
using Random
using Printf  
using Statistics  
import BendersDecomposition: generate_cuts
include("$(dirname(@__DIR__))/example/uflp/data_reader.jl")
include("$(dirname(@__DIR__))/example/uflp/oracle.jl")
include("$(dirname(@__DIR__))/example/uflp/model.jl")
global_logger(ConsoleLogger(stderr, Logging.Debug))

frequency = 500
reuse_dcglp = false
add_benders_cuts_to_master = 2
fraction_of_benders_cuts_to_master = 0.05
threads = 7
split_index_selection_rule = LargestFractional()
time_limit = 14400.0
p = Inf
adjust_t_to_fx = false
lift = false

# load settings
args = parse_commandline()

Random.seed!(args["seed"])
instance = args["instance"]
output_dir = args["output_dir"]

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
problem = read_Simple_data(instance)
dim_x = problem.n_facilities
dim_t = problem.n_customers
c_x = problem.fixed_costs
# c_x = zeros(dim_x)
c_t = ones(dim_t)
data = Data(dim_x, dim_t, problem, c_x, c_t)

# -----------------------------------------------------------------------------
# load parameters
# -----------------------------------------------------------------------------
# Algorithm parameters
benders_param = BendersBnBParam(
    time_limit = time_limit,
    gap_tolerance = 1e-6,
    #disjunctive_root_process = true,
    no_good_cut = false,
    verbose = true
)

dcglp_param = DcglpParam(
    time_limit = 100.0,
    gap_tolerance = 1e-3,
    halt_limit = 3,
    iter_limit = 250,
    verbose = true
)

# Solver parameters
master_solver_param = Dict(
    "solver" => "CPLEX", 
    "CPX_PARAM_EPINT" => 1e-9, 
    "CPX_PARAM_EPRHS" => 1e-9,
    "CPX_PARAM_EPGAP" => 1e-6,
    "CPX_PARAM_BRDIR" => 1
    # "CPXPARAM_MIP_Cuts_Nodecuts" => -1
    # , "CPX_PARAM_REPEATPRESOLVE" => 0
)

dcglp_solver_param = Dict(
    "solver" => "CPLEX", 
    "CPX_PARAM_EPRHS" => 1e-9, 
    "CPX_PARAM_NUMERICALEMPHASIS" => 1, 
    "CPX_PARAM_EPOPT" => 1e-9,
    "CPX_PARAM_THREADS" => threads
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
    norm = LpNorm(p), 
    split_index_selection_rule = split_index_selection_rule,
    disjunctive_cut_append_rule = AllDisjunctiveCuts(), 
    strengthened = true, 
    add_benders_cuts_to_master = add_benders_cuts_to_master,  
    fraction_of_benders_cuts_to_master = fraction_of_benders_cuts_to_master, 
    reuse_dcglp = reuse_dcglp, 
    lift = lift, 
    adjust_t_to_fx = adjust_t_to_fx
) 

set_parameter!(disjunctive_oracle, oracle_param)
update_model!(disjunctive_oracle, data)

# -----------------------------------------------------------------------------
# root node preprocessing
# -----------------------------------------------------------------------------
root_seq_type = BendersSeq
root_param = BendersSeqParam(
    time_limit = 100.0,
    gap_tolerance = 1e-9,
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
user_callback = UserCallback(disjunctive_oracle; params=UserCallbackParam(frequency=frequency))

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





