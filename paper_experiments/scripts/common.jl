using ArgParse
using BendersX
using CPLEX
using CSV
using DataFrames
using Gurobi
using JuMP
using Random
using Statistics

const PAPER_ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW_RESULTS_DIR = joinpath(PAPER_ROOT, "results", "raw")
const PROCESSED_RESULTS_DIR = joinpath(PAPER_ROOT, "results", "processed")

const UFLP_CORRECTNESS_INSTANCES = ["p$(i)" for i in 1:71]
const CFLP_CORRECTNESS_INSTANCES = ["p$(i)" for i in 1:71]
const UFLP_250_INSTANCES = vcat(["ga250c-$(i)" for i in 1:5], ["gs250c-$(i)" for i in 1:5])
const CFLP_200_INSTANCES = vcat(["T200x200_5_$(i)" for i in 1:5], ["T200x200_10_$(i)" for i in 1:5])
const SCFLP_50X100_INSTANCES = vcat(
    ["f50-c100-s128-r5-$(i)" for i in 1:5],
    ["f50-c100-s256-r5-$(i)" for i in 1:5],
    ["f50-c100-s512-r5-$(i)" for i in 1:5],
)
const SCFLP_50X50_S512_INSTANCES = ["f50-c50-s512-r5-$(i)" for i in 1:5]

ensure_dir(path::AbstractString) = (isdir(path) || mkpath(path); path)

const DEFAULT_SOLVER = "gurobi"
const SUPPORTED_SOLVERS = ("gurobi", "cplex")

function normalize_solver_name(solver_name::AbstractString)
    solver = lowercase(strip(solver_name))
    solver in SUPPORTED_SOLVERS || error("Unsupported solver $(solver_name). Expected one of: $(join(SUPPORTED_SOLVERS, ", ")).")
    return solver
end

function lp_optimizer(solver_name::AbstractString; threads::Int = 1, silent::Bool = true)
    solver = normalize_solver_name(solver_name)
    if solver == "gurobi"
        return optimizer_with_attributes(
            Gurobi.Optimizer,
            "Threads" => threads,
            "FeasibilityTol" => 1e-9,
            "OptimalityTol" => 1e-9,
            "NumericFocus" => 1,
            MOI.Silent() => silent,
        )
    elseif solver == "cplex"
        return optimizer_with_attributes(
            CPLEX.Optimizer,
            "CPXPARAM_Threads" => threads,
            "CPX_PARAM_EPRHS" => 1e-9,
            "CPX_PARAM_EPOPT" => 1e-9,
            "CPX_PARAM_NUMERICALEMPHASIS" => 1,
            MOI.Silent() => silent,
        )
    end
end

function mip_optimizer(solver_name::AbstractString; threads::Int = 1, silent::Bool = true, mip_gap::Float64 = 1e-6)
    solver = normalize_solver_name(solver_name)
    if solver == "gurobi"
        return optimizer_with_attributes(
            Gurobi.Optimizer,
            "Threads" => threads,
            "IntFeasTol" => 1e-9,
            "FeasibilityTol" => 1e-9,
            "OptimalityTol" => 1e-9,
            "MIPGap" => mip_gap,
            "NumericFocus" => 1,
            MOI.Silent() => silent,
        )
    elseif solver == "cplex"
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
end

function update_uflp_knapsack_master!(model::Model, data::UFLPData)
    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t[1:data.n_customers] >= -1e6)
    @constraint(model, sum(x) >= 2)
    @objective(model, Min, data.fixed_costs' * x + sum(t))
    return (x = x,), t
end

function generate_stochastic_capacited_facility_location(
    n_facilities::Int,
    n_customers::Int,
    n_scenarios::Int,
    ratio::Int;
    rng::AbstractRNG = Random.default_rng(),
)
    c_x = rand(rng, n_customers)
    c_y = rand(rng, n_customers)

    f_x = rand(rng, n_facilities)
    f_y = rand(rng, n_facilities)

    base_demands = rand(rng, 5:35, n_customers)
    demand_stds = rand(rng, n_customers) .* (0.02 .* base_demands) .+ (0.01 .* base_demands)

    demands = Vector{Vector{Float64}}(undef, n_scenarios)
    demands[1] = Float64.(base_demands)
    for s in 2:n_scenarios
        scenario_demands = zeros(Int, n_customers)
        for j in 1:n_customers
            scenario_demands[j] = max(1, round(Int, base_demands[j] + demand_stds[j] * randn(rng)))
        end
        demands[s] = Float64.(scenario_demands)
    end

    capacities = rand(rng, 10:160, n_facilities)
    fixed_costs = (rand(rng, 100:110, n_facilities) .* sqrt.(capacities)) .+ rand(rng, 0:90, n_facilities)
    fixed_costs = round.(Int, fixed_costs)

    total_demand_max = maximum(sum.(demands))
    total_capacity = sum(capacities)
    capacities = capacities .* ratio .* total_demand_max ./ total_capacity
    capacities = round.(Int, capacities)

    trans_costs = sqrt.((f_x .- c_x') .^ 2 .+ (f_y .- c_y') .^ 2) .* 10

    return SCFLPData(
        n_facilities,
        n_customers,
        n_scenarios,
        Float64.(capacities),
        demands,
        Float64.(fixed_costs),
        Float64.(trans_costs),
    )
end

function load_generated_scflp_instance(instance::AbstractString)
    m = match(r"^f(\d+)-c(\d+)-s(\d+)-r(\d+)-(\d+)$", instance)
    m === nothing && error("Cannot parse SCFLP instance name $(instance). Expected f50-c100-s128-r5-1.")
    n_facilities = parse(Int, m.captures[1])
    n_customers = parse(Int, m.captures[2])
    n_scenarios = parse(Int, m.captures[3])
    ratio = parse(Int, m.captures[4])
    replicate = parse(Int, m.captures[5])
    seed = 1_000_000_000 + n_facilities * 1_000_000 + n_customers * 10_000 + n_scenarios * 100 + ratio * 10 + replicate
    return generate_stochastic_capacited_facility_location(n_facilities, n_customers, n_scenarios, ratio; rng = MersenneTwister(seed))
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
        data = load_generated_scflp_instance(instance)
        return data, "generated_original_scflp"
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
    elseif env_name == "callback" || env_name == "bnb"
        param = BendersBnBParam(; time_limit = time_limit, gap_tolerance = gap_tolerance, verbose = verbose)
        return BendersBnB(master, oracle; param = param)
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
    solver_name::AbstractString,
    summary_file::AbstractString,
    trace_file::AbstractString,
    baseline_obj::Float64 = NaN,
    verbose::Bool = false,
    lp_relaxation::Bool = false,
)
    solver = normalize_solver_name(solver_name)
    master_optimizer = mip_optimizer(solver; threads = solver_threads, silent = !verbose)
    oracle_optimizer = lp_optimizer(solver; threads = solver_threads, silent = !verbose)
    master = build_master_for_oracle(data, oracle_name; optimizer = master_optimizer)
    lp_relaxation && relax_integrality(master.model)
    oracle = build_oracle_for_case(data, master, oracle_name; optimizer = oracle_optimizer)
    env = build_env_for_case(master, oracle, env_name, data; time_limit = time_limit, gap_tolerance = gap_tolerance, verbose = verbose)

    log = DataFrame()
    elapsed = @elapsed begin
        log = solve!(env)
    end

    status = status_string(env.termination_status)
    obj = isfinite(env.obj_value) ? env.obj_value : NaN
    final_lb = hasproperty(log, :LB) ? safe_last(log.LB) : (hasproperty(log, :obj_bound) ? safe_last(log.obj_bound) : NaN)
    final_ub = hasproperty(log, :UB) ? safe_last(log.UB) : (hasproperty(log, :obj_val) ? safe_last(log.obj_val) : NaN)
    final_gap = hasproperty(log, :gap) ? safe_last(log.gap) : (hasproperty(log, :rel_gap) ? safe_last(log.rel_gap) : NaN)
    master_time = hasproperty(log, :master_time) ? safe_sum(log.master_time) : NaN
    oracle_time = hasproperty(log, :oracle_time) ? safe_sum(log.oracle_time) : NaN
    total_iter_time = hasproperty(log, :total_time) ? safe_sum(log.total_time) : (hasproperty(log, :time) ? safe_last(log.time) : elapsed)
    oracle_share = isfinite(oracle_time) && total_iter_time > 0 ? oracle_time / total_iter_time : NaN
    iterations = hasproperty(log, :node_count) ? safe_last(log.node_count) : nrow(log)
    root_node_time = hasproperty(log, :root_node_time) ? safe_last(log.root_node_time) : NaN
    node_count = hasproperty(log, :node_count) ? safe_last(log.node_count) : missing
    n_lazy_cuts = hasproperty(log, :n_lazy_cuts) ? safe_last(log.n_lazy_cuts) : missing
    n_user_cuts = hasproperty(log, :n_user_cuts) ? safe_last(log.n_user_cuts) : missing
    dims = size_metadata(data)
    error = isfinite(baseline_obj) && isfinite(obj) ? abs(obj - baseline_obj) : NaN

    metadata = (
        experiment = experiment,
        config_type = lp_relaxation ? "benders_lp_relax" : "benders",
        problem = problem_label(problem),
        instance = instance,
        data_source = data_source,
        n_facilities = dims.n_facilities,
        n_customers = dims.n_customers,
        n_scenarios = dims.n_scenarios,
        env = env_name,
        oracle = oracle_name,
        repeat = repeat,
        solver = solver,
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
        iterations = iterations,
        root_node_time = root_node_time,
        node_count = node_count,
        n_lazy_cuts = n_lazy_cuts,
        n_user_cuts = n_user_cuts,
        baseline_obj = baseline_obj,
        obj_error = error,
    ))
    append_csv_row!(summary_file, row)
    return row
end

function solve_extensive_mip(problem::AbstractString, data; time_limit::Float64, solver_threads::Int, solver_name::AbstractString, verbose::Bool = false)
    solver = normalize_solver_name(solver_name)
    model = Model(mip_optimizer(solver; threads = solver_threads, silent = !verbose))
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
    solver_name::AbstractString,
    summary_file::AbstractString,
    verbose::Bool = false,
)
    solver = normalize_solver_name(solver_name)
    result = solve_extensive_mip(problem, data; time_limit = time_limit, solver_threads = solver_threads, solver_name = solver, verbose = verbose)
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
        solver = solver,
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
