using BendersX
using JuMP
using CPLEX
using Printf
using Random

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))
include(normpath(joinpath(@__DIR__, "..", "script_utils.jl")))


options, _ = parse_script_args(ARGS)

instance = get_string_option(options, "instance", "T100x100_5_2")
time_limit = get_float_option(options, "time_limit", 60.0)
threads = get_int_option(options, "threads", 1)
values_arg = get_string_option(options, "values", "both")
silent = get_bool_option(options, "silent", false)
seed = get_int_option(options, "seed", 1)
output_csv = get_string_option(options, "output_csv", "")

Random.seed!(seed)

function build_cflp_milp(data::CFLPData)
    model = Model(mip_optimizer)
    I, J = data.n_facilities, data.n_customers

    @variable(model, x[1:I], Bin)
    @variable(model, y[1:I, 1:J] >= 0)
    @variable(model, t)

    cost_demands = data.costs .* data.demands'
    @objective(model, Min, data.fixed_costs' * x + t)

    @constraint(model, recourse_cost, t >= sum(cost_demands .* y))
    @constraint(model, demand[j in 1:J], sum(y[:, j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.demands .* y[i, :]) <= data.capacities[i] * x[i])
    @constraint(model, capacity_total, sum(data.capacities[i] * x[i] for i in 1:I) >= sum(data.demands))

    return model, x
end

function configure_solver!(model::Model, threads::Int, time_limit::Float64, silent::Bool)
    set_optimizer_attribute(model, "CPXPARAM_Threads", threads)
    set_time_limit_sec(model, time_limit)
    set_optimizer_attribute(model, MOI.Silent(), silent)
    return nothing
end

function objective_or_nan(model::Model)
    return has_values(model) ? objective_value(model) : NaN
end

function bound_or_nan(model::Model)
    try
        return objective_bound(model)
    catch
        return NaN
    end
end

function values_to_scan(values::String)
    if values == "both"
        return [0, 1]
    elseif values == "zero"
        return [0]
    elseif values == "one"
        return [1]
    end
    throw(ArgumentError("Unsupported `values` option `$(values)`. Use both, zero, or one."))
end

function scan_fixings!(model::Model, x, fix_values::Vector{Int}, max_variables::Int)
    n = length(x)
    n_scan = max_variables > 0 ? min(max_variables, n) : n
    rows = NamedTuple[]

    for i in 1:n_scan
        for v in fix_values
            fix(x[i], v; force = true)
            optimize!(model)

            status = termination_status(model)
            pstatus = JuMP.primal_status(model)

            push!(
                rows,
                (
                    variable_index = i,
                    fixed_value = v,
                    termination_status = string(status),
                    primal_status = string(pstatus),
                    solve_time = solve_time(model),
                    objective_value = objective_or_nan(model),
                    objective_bound = bound_or_nan(model),
                    is_infeasible = status in (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED),
                ),
            )

            unfix(x[i])
        end
    end

    return rows
end

function write_results_csv(path::String, rows::Vector{<:NamedTuple})
    headers = (
        :variable_index,
        :fixed_value,
        :termination_status,
        :primal_status,
        :solve_time,
        :objective_value,
        :objective_bound,
        :is_infeasible,
    )
    open(path, "w") do io
        println(io, join(string.(headers), ","))
        for row in rows
            values = [string(getproperty(row, key)) for key in headers]
            println(io, join(values, ","))
        end
    end
end

function print_infeasible_rows(rows::Vector{<:NamedTuple})
    for row in rows
        @printf(
            "  x[%d] = %d -> status=%s, solve_time=%.4f\n",
            row.variable_index,
            row.fixed_value,
            row.termination_status,
            row.solve_time,
        )
    end
end

data = read_cfl_file(instance)
model, x = build_cflp_milp(data)
configure_solver!(model, threads, time_limit, silent)

@printf("CFLP MILP fix scan\n")
@printf("  instance        : %s\n", instance)
@printf("  n_facilities    : %d\n", data.n_facilities)
@printf("  n_customers     : %d\n", data.n_customers)
@printf("  time_limit      : %.2f\n", time_limit)
@printf("  threads         : %d\n", threads)

optimize!(model)
base_status = termination_status(model)
@printf("  baseline_status : %s\n", string(base_status))
if has_values(model)
    @printf("  baseline_obj    : %.6f\n", objective_value(model))
end

fix_values = values_to_scan(values_arg)
results = scan_fixings!(model, x, fix_values, 0)
infeasible_rows = filter(row -> row.is_infeasible, results)

@printf("  scanned_fixes   : %d\n", length(results))
@printf("  infeasible_fixes: %d\n", length(infeasible_rows))

if !isempty(infeasible_rows)
    println("\nInfeasible variable fixings:")
    print_infeasible_rows(infeasible_rows)
end

if !isempty(output_csv)
    write_results_csv(output_csv, results)
end
