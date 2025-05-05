using JuMP, DataFrames, Logging, CSV
using BendersDecomposition
using Printf  
using Statistics  
import BendersDecomposition: generate_cuts
include("$(dirname(@__DIR__))/example/cflp/data_reader.jl")
include("$(dirname(@__DIR__))/example/cflp/oracle.jl")
include("$(dirname(@__DIR__))/example/cflp/model.jl")

# load settings
args = parse_commandline()

instance = args["instance"]
output_dir = args["output_dir"]

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
master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 4)
typical_oracal_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)
# -----------------------------------------------------------------------------
# MIP model
# -----------------------------------------------------------------------------
master = Master(data; solver_param = master_solver_param)
update_model!(master, data)

oracle = CFLKnapsackOracle(data; solver_param = typical_oracal_solver_param)
update_model!(oracle, data)

env = BendersSeqInOut(data, master, oracle; param = benders_inout_param)
relax_integrality(env.master.model)

log = solve!(env)

@info termination_status(env.master.model)
@info "Solve time: $(solve_time(env.master.model))"
@info "Objective value: $(objective_value(env.master.model))"
@info "objective bound: $(objective_bound(env.master.model))"