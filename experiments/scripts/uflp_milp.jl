using JuMP, DataFrames, Logging, CSV
using BendersX
using ArgParse
using Random
using Printf
using Statistics
include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--instance"
            help = "Instance name"
            arg_type = String
            default = "ga250a-1"
        "--seed"
            help = "Random seed"
            arg_type = Int
            default = 1
        "--output_dir"
            help = "Output directory"
            arg_type = String
            default = "output"
    end
    return parse_args(s)
end

function create_mip_model(data::UFLPData)
    model = Model(mip_optimizer)
    I, J = data.n_facilities, data.n_customers
    @variable(model, x[1:I], Bin)
    @variable(model, y[1:I, 1:J] >= 0)
    @variable(model, t[1:J] >= 0)

    cost_demands = data.costs .* data.demands'
    @objective(model, Min, data.fixed_costs' * x + sum(t))

    @constraint(model, obj[j in 1:J], t[j] >= sum(cost_demands[:, j] .* y[:, j]))
    @constraint(model, demand[j in 1:J], sum(y[:, j]) == 1)
    @constraint(model, facility_open, y .<= x)
    return model
end

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
data = read_Simple_data(instance)

# -----------------------------------------------------------------------------
# MIP model
# -----------------------------------------------------------------------------
mip_model = create_mip_model(data)
set_optimizer_attribute(mip_model, "CPXPARAM_Threads", threads)
set_optimizer_attribute(mip_model, "CPX_PARAM_BRDIR", 1)
set_time_limit_sec(mip_model, time_limit)
set_optimizer_attribute(mip_model, MOI.Silent(), false)
optimize!(mip_model)

@info termination_status(mip_model)
@info "Node count: $(node_count(mip_model))"
@info "Elapsed time: $(solve_time(mip_model))"
@info "Objective value: $(objective_value(mip_model))"
@info "Objective bound: $(objective_bound(mip_model))"
@info "Relative gap: $(relative_gap(mip_model))"
