using JuMP, DataFrames, Logging, CSV
using BendersX
using ArgParse
using Printf
using Statistics
include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--instance"
            help = "Instance name"
            arg_type = String
            default = "T100x100_5_1"
        "--output_dir"
            help = "Output directory"
            arg_type = String
            default = "output"
    end
    return parse_args(s)
end

function create_mip_model(data::CFLPData)
    model = Model(mip_optimizer)
    I, J = data.n_facilities, data.n_customers
    @variable(model, x[1:I], Bin)
    @variable(model, y[1:I, 1:J] >= 0)
    @variable(model, t)

    cost_demands = data.costs .* data.demands'
    @objective(model, Min, data.fixed_costs' * x + t)

    @constraint(model, t >= sum(cost_demands .* y))
    @constraint(model, demand[j in 1:J], sum(y[:, j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.demands .* y[i, :]) <= data.capacities[i] * x[i])
    @constraint(model, capacity_total, sum(data.capacities[i] * x[i] for i in 1:I) >= sum(data.demands))
    return model
end

# load settings
args = parse_commandline()

instance = args["instance"]
output_dir = args["output_dir"]

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
data = read_cfl_file(instance)

# -----------------------------------------------------------------------------
# MIP model
# -----------------------------------------------------------------------------
mip_model = create_mip_model(data)
set_optimizer_attribute(mip_model, "CPXPARAM_Threads", 7)
set_time_limit_sec(mip_model, 14400.0)
set_optimizer_attribute(mip_model, MOI.Silent(), false)
optimize!(mip_model)

@info termination_status(mip_model)
@info "Solve time: $(solve_time(mip_model))"
@info "Objective value: $(objective_value(mip_model))"
@info "Objective bound: $(objective_bound(mip_model))"
