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
    disjunctive_root_process = false,
    verbose = true
)

dcglp_param = DcglpParam(
    time_limit = 200.0,
    gap_tolerance = 1e-3,
    halt_limit = 3,
    iter_limit = 15,
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

dcglp_solver_param = Dict(
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
typical_oracles = [
    SeparableOracle(data, ClassicalOracle(), data.problem.num_scenarios; solver_param = typical_oracle_solver_param),
    SeparableOracle(data, ClassicalOracle(), data.problem.num_scenarios; solver_param = typical_oracle_solver_param)
]

for k in 1:2, j in 1:typical_oracles[k].N
    update_model!(typical_oracles[k].oracles[j], data, j)
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
    split_index_selection_rule = LargestFractional(),
    disjunctive_cut_append_rule = AllDisjunctiveCuts(), 
    strengthened = true, 
    add_benders_cuts_to_master = true, 
    fraction_of_benders_cuts_to_master = 1.0, 
    reuse_dcglp = true,
    lift = true,
    adjust_t_to_fx = false
)

set_parameter!(disjunctive_oracle, oracle_param)
update_model!(disjunctive_oracle, data)

# -----------------------------------------------------------------------------
# root node preprocessing
# -----------------------------------------------------------------------------
# root_seq_type = BendersSeqInOut
# root_param = BendersSeqInOutParam(
#     time_limit = 1800.0,
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

lazy_oracle = SeparableOracle(data, ClassicalOracle(), data.problem.num_scenarios; solver_param = typical_oracle_solver_param)
for j=1:lazy_oracle.N
    update_model!(lazy_oracle.oracles[j], data, j)
end

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






