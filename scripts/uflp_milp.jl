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

threads = 7
time_limit = 14400.0

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
c_t = ones(dim_t)
data = Data(dim_x, dim_t, problem, c_x, c_t)

# -----------------------------------------------------------------------------
# load parameters
# -----------------------------------------------------------------------------
mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, 
"CPXPARAM_Threads" => threads,  "CPX_PARAM_BRDIR" => 1)

# -----------------------------------------------------------------------------
# MIP model
# -----------------------------------------------------------------------------
mip = Mip(data)
assign_attributes!(mip.model, mip_solver_param)
set_time_limit_sec(mip.model, time_limit)
update_model!(mip, data)
set_optimizer_attribute(mip.model, MOI.Silent(), false)
optimize!(mip.model)

@info termination_status(mip.model)
@info "Node count: $(node_count(mip.model))"
@info "Elapsed time: $(solve_time(mip.model))"
@info "Objective value: $(objective_value(mip.model))"
@info "Objective bound: $(objective_bound(mip.model))"
@info "Relative gap: $(relative_gap(mip.model))"
