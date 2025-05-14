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
mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 32)

# -----------------------------------------------------------------------------
# MIP model
# -----------------------------------------------------------------------------
mip = Mip(data)
assign_attributes!(mip.model, mip_solver_param)
set_time_limit_sec(mip.model, 3600.0)
update_model!(mip, data)
set_optimizer_attribute(mip.model, MOI.Silent(), false)
optimize!(mip.model)

@info termination_status(mip.model)
@info "Solve time: $(solve_time(mip.model))"
@info "Objective value: $(objective_value(mip.model))"
@info "objective bound: $(objective_bound(mip.model))"