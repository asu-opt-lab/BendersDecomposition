using BendersDecomposition
using JuMP
using CPLEX
using Gurobi
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
# appy parameters
# -----------------------------------------------------------------------------
mip_solver_param = Dict("solver" => "CPLEX", "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6)

assign_attributes!(model, mip_solver_param)
set_time_limit_sec(model, 3600.0)
set_optimizer_attribute(model, MOI.Silent(), false)
optimize!(model)

# -----------------------------------------------------------------------------
# solve
# -----------------------------------------------------------------------------
println("termination status: $(termination_status(model))")
println("objective bound: $(objective_bound(model))")
println("Objective value: $(objective_value(model))")