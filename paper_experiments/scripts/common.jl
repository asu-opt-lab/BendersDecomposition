using ArgParse
using BendersX
using CPLEX
using CSV
using DataFrames
using JuMP
using Random
using Statistics

const PAPER_ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW_RESULTS_DIR = joinpath(PAPER_ROOT, "results", "raw")
const PROCESSED_RESULTS_DIR = joinpath(PAPER_ROOT, "results", "processed")

const UFLP_CORRECTNESS_INSTANCES = ["p$(i)" for i in 1:71]
const CFLP_CORRECTNESS_INSTANCES = ["p$(i)" for i in 1:71]
const UFLP_250_INSTANCES = vcat(["ga250a-$(i)" for i in 1:5], ["ga250b-$(i)" for i in 1:5])
const CFLP_200_INSTANCES = vcat(["T200x200_5_$(i)" for i in 1:5], ["T200x200_10_$(i)" for i in 1:5])
const SCFLP_100X200_INSTANCES = ["synthetic-f100-c200-s256-$(i)" for i in 1:10]

ensure_dir(path::AbstractString) = (isdir(path) || mkpath(path); path)

function cplex_lp_optimizer(; threads::Int = 1, silent::Bool = true)
    return optimizer_with_attributes(
        CPLEX.Optimizer,
        "CPXPARAM_Threads" => threads,
        "CPX_PARAM_EPRHS" => 1e-9,
        "CPX_PARAM_EPOPT" => 1e-9,
        "CPX_PARAM_NUMERICALEMPHASIS" => 1,
        MOI.Silent() => silent,
    )
end

function cplex_mip_optimizer(; threads::Int = 1, silent::Bool = true, mip_gap::Float64 = 1e-6)
    return optimizer_with_attributes(
        CPLEX.Optimizer,
        "CPXPARAM_Threads" => threads,
        "CPX_PARAM_EPINT" => 1e-9,
        "CPX_PARAM_EPRHS" => 1e-9,
        "CPX_PARAM_EPGAP" => mip_gap,
        "CPX_PARAM_EPOPT" => 1e-9,
        "CPX_PARAM_NUMERICALEMPHASIS" => 1,
        MOI.Silent() => silent,
    )
end

function update_uflp_knapsack_master!(model::Model, data::UFLPData)
    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t[1:data.n_customers] >= -1e6)
    @constraint(model, sum(x) >= 2)
    @objective(model, Min, data.fixed_costs' * x + sum(t))
    return (x = x,), t
end

function synthetic_scflp_data(
    instance_id::Int;
    n_facilities::Int = 100,
    n_customers::Int = 200,
    n_scenarios::Int = 256,
    seed_base::Int = 430_000,
)
    rng = MersenneTwister(seed_base + instance_id)
    facility_xy = rand(rng, n_facilities, 2)
    customer_xy = rand(rng, n_customers, 2)

    costs = zeros(Float64, n_facilities, n_customers)
    for i in 1:n_facilities, j in 1:n_customers
        dx = facility_xy[i, 1] - customer_xy[j, 1]
        dy = facility_xy[i, 2] - customer_xy[j, 2]
        costs[i, j] = 10.0 + 100.0 * sqrt(dx * dx + dy * dy) + rand(rng)
    end

    fixed_costs = Float64.(rand(rng, 800:2600, n_facilities))
    capacities = Float64.(rand(rng, 70:110, n_facilities))
    demands = [Float64.(rand(rng, 10:24, n_customers)) for _ in 1:n_scenarios]

    max_demand = maximum(sum.(demands))
    if sum(capacities) < 1.35 * max_demand
        capacities .*= 1.35 * max_demand / sum(capacities)
    end

    return SCFLPData(n_facilities, n_customers, n_scenarios, capacities, demands, fixed_costs, costs)
end

function parse_synthetic_scflp_index(instance::AbstractString)
    m = match(r"-(\d+)$", instance)
    m === nothing && error("Cannot parse synthetic SCFLP instance index from $(instance)")
    return parse(Int, m.captures[1])
end

function load_paper_instance(problem::AbstractString, instance::AbstractString)
    if problem == "uflp_correctness"
        data = read_uflp_benchmark_data(instance)
        return data, "uflp_locssall"
    elseif problem == "cflp_correctness"
        data = read_cflp_benchmark_data(instance)
        return data, "cflp_locssall"
    elseif problem == "uflp"
        data = read_Simple_data(instance)
        return data, "uflp_allkoerkelghosh"
    elseif problem == "cflp"
        data = read_cfl_file(instance)
        return data, "cflp_output"
    elseif problem == "scflp"
        data = synthetic_scflp_data(parse_synthetic_scflp_index(instance))
        return data, "synthetic_scflp"
    else
        error("Unknown problem $(problem)")
    end
end

problem_label(problem::AbstractString) = replace(problem, "_correctness" => "")

function size_metadata(data)
    if data isa UFLPData
        return (n_facilities = data.n_facilities, n_customers = data.n_customers, n_scenarios = 0)
    elseif data isa CFLPData
        return (n_facilities = data.n_facilities, n_customers = data.n_customers, n_scenarios = 0)
    elseif data isa SCFLPData
        return (n_facilities = data.n_facilities, n_customers = data.n_customers, n_scenarios = data.n_scenarios)
    else
        return (n_facilities = 0, n_customers = 0, n_scenarios = 0)
    end
end

function git_commit()
    try
        return readchomp(`git rev-parse --short HEAD`)
    catch
        return "unknown"
    end
end

function append_csv_row!(path::AbstractString, row::NamedTuple)
    ensure_dir(dirname(path))
    incoming = DataFrame([row])
    if !isfile(path)
        CSV.write(path, incoming)
        return path
    end

    existing = DataFrame(CSV.File(path))
    column_order = vcat(names(existing), setdiff(names(incoming), names(existing)))
    for name in setdiff(column_order, names(existing))
        existing[!, name] = fill(missing, nrow(existing))
    end
    for name in setdiff(column_order, names(incoming))
        incoming[!, name] = [missing]
    end
    CSV.write(path, vcat(existing[:, column_order], incoming[:, column_order]; cols = :setequal))
    return path
end

function append_trace!(path::AbstractString, log::DataFrame, metadata::NamedTuple)
    ensure_dir(dirname(path))
    df = copy(log)
    for (k, v) in pairs(metadata)
        df[!, Symbol(k)] = fill(v, nrow(df))
    end
    CSV.write(path, df; append = isfile(path), writeheader = !isfile(path))
    return path
end

safe_last(xs, default = NaN) = isempty(xs) ? default : xs[end]
safe_sum(xs) = isempty(xs) ? 0.0 : sum(skipmissing(xs))

function status_string(status)
    return string(status)
end

is_optimal_status(status::AbstractString) = occursin("Optimal", status)

function build_master_for_oracle(data, oracle_name::AbstractString; optimizer)
    builder = (data isa UFLPData && oracle_name == "ufl_knapsack") ? update_uflp_knapsack_master! : update_master_model!
    return Master(data; model = builder, optimizer = optimizer)
end

function build_oracle_for_case(data, master, oracle_name::AbstractString; optimizer)
    if data isa SCFLPData
        if oracle_name == "classical"
            return SeparableOracle(data, master, ClassicalOracle(), data.n_scenarios; model = update_sub_model!, optimizer = optimizer)
        elseif oracle_name == "unified"
            return SeparableOracle(data, master, UnifiedOracle(), data.n_scenarios; model = update_sub_model!, sub_oracle_param = UnifiedOracleParam(), optimizer = optimizer)
        elseif oracle_name == "pareto"
            param = ParetoOracleParam(fill(1.0, data.n_facilities))
            return SeparableOracle(data, master, ParetoOracle(), data.n_scenarios; model = update_sub_model!, sub_oracle_param = param, optimizer = optimizer)
        elseif oracle_name == "cfl_knapsack"
            return SeparableOracle(data, master, CFLKnapsackOracle(), data.n_scenarios; model = update_sub_model!, optimizer = optimizer)
        else
            error("Unsupported SCFLP oracle $(oracle_name)")
        end
    end

    if oracle_name == "classical"
        return ClassicalOracle(data, master; model = update_sub_model!, optimizer = optimizer)
    elseif oracle_name == "unified"
        return UnifiedOracle(data, master; model = update_sub_model!, optimizer = optimizer)
    elseif oracle_name == "pareto"
        param = ParetoOracleParam(fill(1.0, master.dim_x))
        return ParetoOracle(data, master, param; model = update_sub_model!, optimizer = optimizer)
    elseif oracle_name == "ufl_knapsack"
        oracle = UFLKnapsackOracle(data)
        set_parameter!(oracle, "add_only_violated_cuts", true)
        return oracle
    elseif oracle_name == "cfl_knapsack"
        return CFLKnapsackOracle(data, master; model = update_sub_model!, optimizer = optimizer)
    else
        error("Unsupported oracle $(oracle_name)")
    end
end

function build_env_for_case(master, oracle, env_name::AbstractString, data; time_limit::Float64, gap_tolerance::Float64, verbose::Bool)
    if env_name == "seq"
        param = BendersSeqParam(; time_limit = time_limit, gap_tolerance = gap_tolerance, verbose = verbose)
        return BendersSeq(master, oracle; param = param)
    elseif env_name == "seq_inout"
        param = BendersSeqInOutParam(;
            time_limit = time_limit,
            gap_tolerance = gap_tolerance,
            verbose = verbose,
            stabilizing_x = ones(master.dim_x),
            α = 0.9,
            λ = 0.1,
        )
        return BendersSeqInOut(master, oracle; param = param)
    else
        error("Unsupported environment $(env_name)")
    end
end

function run_benders_case(;
    experiment::AbstractString,
    problem::AbstractString,
    instance::AbstractString,
    data,
    data_source::AbstractString,
    env_name::AbstractString,
    oracle_name::AbstractString,
    repeat::Int,
    time_limit::Float64,
    gap_tolerance::Float64,
    solver_threads::Int,
    summary_file::AbstractString,
    trace_file::AbstractString,
    baseline_obj::Float64 = NaN,
    verbose::Bool = false,
)
    master_optimizer = cplex_mip_optimizer(; threads = solver_threads, silent = !verbose)
    oracle_optimizer = cplex_lp_optimizer(; threads = solver_threads, silent = !verbose)
    master = build_master_for_oracle(data, oracle_name; optimizer = master_optimizer)
    oracle = build_oracle_for_case(data, master, oracle_name; optimizer = oracle_optimizer)
    env = build_env_for_case(master, oracle, env_name, data; time_limit = time_limit, gap_tolerance = gap_tolerance, verbose = verbose)

    log = DataFrame()
    elapsed = @elapsed begin
        log = solve!(env)
    end

    status = status_string(env.termination_status)
    obj = isfinite(env.obj_value) ? env.obj_value : NaN
    final_lb = hasproperty(log, :LB) ? safe_last(log.LB) : NaN
    final_ub = hasproperty(log, :UB) ? safe_last(log.UB) : NaN
    final_gap = hasproperty(log, :gap) ? safe_last(log.gap) : NaN
    master_time = hasproperty(log, :master_time) ? safe_sum(log.master_time) : NaN
    oracle_time = hasproperty(log, :oracle_time) ? safe_sum(log.oracle_time) : NaN
    total_iter_time = hasproperty(log, :total_time) ? safe_sum(log.total_time) : elapsed
    oracle_share = total_iter_time > 0 ? oracle_time / total_iter_time : NaN
    dims = size_metadata(data)
    error = isfinite(baseline_obj) && isfinite(obj) ? abs(obj - baseline_obj) : NaN

    metadata = (
        experiment = experiment,
        config_type = "benders",
        problem = problem_label(problem),
        instance = instance,
        data_source = data_source,
        n_facilities = dims.n_facilities,
        n_customers = dims.n_customers,
        n_scenarios = dims.n_scenarios,
        env = env_name,
        oracle = oracle_name,
        repeat = repeat,
        julia_threads = Threads.nthreads(),
        solver_threads = solver_threads,
        time_limit = time_limit,
        gap_tolerance = gap_tolerance,
        git_commit = git_commit(),
    )
    append_trace!(trace_file, log, metadata)

    row = merge(metadata, (
        status = status,
        solved = is_optimal_status(status),
        obj_value = obj,
        LB = final_lb,
        UB = final_ub,
        final_gap = final_gap,
        total_time = elapsed,
        iteration_time = total_iter_time,
        master_time = master_time,
        oracle_time = oracle_time,
        oracle_share = oracle_share,
        iterations = nrow(log),
        baseline_obj = baseline_obj,
        obj_error = error,
    ))
    append_csv_row!(summary_file, row)
    return row
end

function solve_extensive_mip(problem::AbstractString, data; time_limit::Float64, solver_threads::Int, verbose::Bool = false)
    model = Model(cplex_mip_optimizer(; threads = solver_threads, silent = !verbose))
    set_time_limit_sec(model, time_limit)

    if data isa UFLPData
        I, J = data.n_facilities, data.n_customers
        @variable(model, x[1:I], Bin)
        @variable(model, y[1:I, 1:J] >= 0)
        cost_demands = data.costs .* data.demands'
        @objective(model, Min, data.fixed_costs' * x + sum(cost_demands .* y))
        @constraint(model, demand[j in 1:J], sum(y[:, j]) == 1)
        @constraint(model, facility_open[i in 1:I, j in 1:J], y[i, j] <= x[i])
    elseif data isa CFLPData
        I, J = data.n_facilities, data.n_customers
        @variable(model, x[1:I], Bin)
        @variable(model, y[1:I, 1:J] >= 0)
        cost_demands = data.costs .* data.demands'
        @objective(model, Min, data.fixed_costs' * x + sum(cost_demands .* y))
        @constraint(model, demand[j in 1:J], sum(y[:, j]) == 1)
        @constraint(model, facility_open[i in 1:I, j in 1:J], y[i, j] <= x[i])
        @constraint(model, capacity[i in 1:I], sum(data.demands[j] * y[i, j] for j in 1:J) <= data.capacities[i] * x[i])
    else
        error("Extensive MIP baseline is only implemented for UFLP and CFLP")
    end

    elapsed = @elapsed optimize!(model)
    status = string(termination_status(model))
    obj = has_values(model) ? objective_value(model) : NaN
    bound = try
        objective_bound(model)
    catch
        NaN
    end
    gap = try
        has_values(model) ? relative_gap(model) : NaN
    catch
        NaN
    end
    nodes = try
        node_count(model)
    catch
        missing
    end
    return (status = status, obj_value = obj, obj_bound = bound, rel_gap = gap, total_time = elapsed, node_count = nodes)
end

function log_baseline_case(;
    experiment::AbstractString,
    problem::AbstractString,
    instance::AbstractString,
    data,
    data_source::AbstractString,
    repeat::Int,
    time_limit::Float64,
    solver_threads::Int,
    summary_file::AbstractString,
    verbose::Bool = false,
)
    result = solve_extensive_mip(problem, data; time_limit = time_limit, solver_threads = solver_threads, verbose = verbose)
    dims = size_metadata(data)
    row = (
        experiment = experiment,
        config_type = "baseline",
        problem = problem_label(problem),
        instance = instance,
        data_source = data_source,
        n_facilities = dims.n_facilities,
        n_customers = dims.n_customers,
        n_scenarios = dims.n_scenarios,
        env = "extensive_form",
        oracle = "mip",
        repeat = repeat,
        julia_threads = Threads.nthreads(),
        solver_threads = solver_threads,
        time_limit = time_limit,
        gap_tolerance = NaN,
        git_commit = git_commit(),
        status = result.status,
        solved = occursin("OPTIMAL", result.status),
        obj_value = result.obj_value,
        LB = result.obj_bound,
        UB = result.obj_value,
        final_gap = result.rel_gap,
        total_time = result.total_time,
        iteration_time = NaN,
        master_time = NaN,
        oracle_time = NaN,
        oracle_share = NaN,
        iterations = missing,
        baseline_obj = result.obj_value,
        obj_error = 0.0,
        node_count = result.node_count,
    )
    append_csv_row!(summary_file, row)
    return row
end

function default_output_files(experiment_name::AbstractString, output_dir::AbstractString)
    raw_dir = ensure_dir(output_dir)
    return (
        summary = joinpath(raw_dir, "$(experiment_name)_summary.csv"),
        trace = joinpath(raw_dir, "$(experiment_name)_trace.csv"),
    )
end
