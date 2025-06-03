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
problem = read_cfl_file(instance)

dim_x = problem.n_facilities
dim_t = 1
c_x = problem.fixed_costs
c_t = [1]
data = Data(dim_x, dim_t, problem, c_x, c_t)

# -----------------------------------------------------------------------------
# load parameters
# -----------------------------------------------------------------------------
# mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6)
mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPGAP" => 1e-6, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9)
# mip_solver_param = Dict("solver" => "Gurobi")
# -----------------------------------------------------------------------------
# MIP model
# -----------------------------------------------------------------------------
mip = Mip(data)
assign_attributes!(mip.model, mip_solver_param)
set_time_limit_sec(mip.model, 14400.0)
update_model!(mip, data)
set_optimizer_attribute(mip.model, MOI.Silent(), false)
optimize!(mip.model)

@info termination_status(mip.model)
@info "Solve time: $(solve_time(mip.model))"
@info "Objective value: $(objective_value(mip.model))"
@info "objective bound: $(objective_bound(mip.model))"