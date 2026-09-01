using ArgParse
using BendersX
using CSV
using DataFrames
using Dates
using JuMP
using Random

include(joinpath(@__DIR__, "solver_composition_factories.jl"))

function parse_commandline()
    settings = ArgParseSettings(
        description = "Run one UFLP instance with one master/subproblem solver composition.",
    )
    @add_arg_table! settings begin
        "--instance"
            help = "Koerkel-Ghosh UFLP instance, for example ga250a-1"
            arg_type = String
            default = "ga250a-1"
        "--master_solver"
            help = "Master MILP solver: cplex or gurobi"
            arg_type = String
            required = true
        "--sub_solver"
            help = "Subproblem LP solver: cplex or gurobi"
            arg_type = String
            required = true
        "--sub_method"
            help = "Subproblem LP method: dual_simplex or default"
            arg_type = String
            default = "dual_simplex"
        "--seed"
            help = "Julia and solver random seed"
            arg_type = Int
            default = 1
        "--threads"
            help = "Threads per solver model; use 1 for the attribution experiment"
            arg_type = Int
            default = 1
        "--time_limit"
            help = "Benders branch-and-bound time limit in seconds"
            arg_type = Float64
            default = 3600.0
        "--gap_tolerance"
            help = "Relative optimality-gap tolerance"
            arg_type = Float64
            default = 1e-6
        "--numeric_tolerance"
            help = "Common feasibility, optimality, and integrality tolerance"
            arg_type = Float64
            default = 1e-8
        "--output_dir"
            help = "Experiment output directory"
            arg_type = String
            default = "output"
        "--skip_warmup"
            help = "Skip the discarded p1 warm-up solve"
            action = :store_true
        "--overwrite"
            help = "Overwrite an existing result CSV for this configuration"
            action = :store_true
    end
    return parse_args(settings)
end

function optimizer_version(model::Model)
    try
        return string(MOI.get(JuMP.backend(model), MOI.SolverVersion()))
    catch
        return "unknown"
    end
end

function build_environment(
    data::UFLPData,
    master_solver::Symbol,
    sub_solver::Symbol,
    sub_method::Symbol;
    threads::Int,
    seed::Int,
    time_limit::Float64,
    gap_tolerance::Float64,
    numeric_tolerance::Float64,
)
    master_optimizer = composition_optimizer(
        master_solver,
        :master;
        threads = threads,
        seed = seed,
        sub_method = sub_method,
        tolerance = numeric_tolerance,
    )
    subproblem_optimizer = composition_optimizer(
        sub_solver,
        :subproblem;
        threads = threads,
        seed = seed,
        sub_method = sub_method,
        tolerance = numeric_tolerance,
    )

    master = Master(
        data;
        model = update_master_model!,
        optimizer = master_optimizer,
    )
    oracle = ClassicalOracle(
        data,
        master;
        model = update_sub_model!,
        optimizer = subproblem_optimizer,
    )
    parameters = BendersBnBParam(
        time_limit = time_limit,
        gap_tolerance = gap_tolerance,
        verbose = false,
    )
    return BendersBnB(master, oracle; param = parameters), oracle
end

function warmup!(
    master_solver::Symbol,
    sub_solver::Symbol,
    sub_method::Symbol;
    threads::Int,
    seed::Int,
    gap_tolerance::Float64,
    numeric_tolerance::Float64,
)
    @info "Warming solver-composition code path; result will be discarded" master_solver sub_solver
    data = read_uflp_benchmark_data("p1")
    env, _ = build_environment(
        data,
        master_solver,
        sub_solver,
        sub_method;
        threads = threads,
        seed = seed,
        time_limit = 60.0,
        gap_tolerance = gap_tolerance,
        numeric_tolerance = numeric_tolerance,
    )
    solve!(env)
    env.termination_status == Optimal() || error(
        "Warm-up solve failed with status $(env.termination_status)",
    )
    return nothing
end

function output_path(
    output_dir::String,
    instance::String,
    master_solver::Symbol,
    sub_solver::Symbol,
    sub_method::Symbol,
    seed::Int,
)
    filename = join(
        (
            instance,
            "master-$(master_solver)",
            "sub-$(sub_solver)",
            "method-$(sub_method)",
            "seed-$(seed)",
        ),
        "__",
    ) * ".csv"
    return joinpath(abspath(output_dir), "results_csv", filename)
end

function write_result(path::String, row; overwrite::Bool)
    mkpath(dirname(path))
    isfile(path) && !overwrite && error(
        "Result file already exists: $path. Pass --overwrite to replace it.",
    )
    CSV.write(path, DataFrame([row]))
    return nothing
end

function main()
    args = parse_commandline()
    master_solver = parse_solver_name(args["master_solver"])
    sub_solver = parse_solver_name(args["sub_solver"])
    sub_method = parse_sub_method(args["sub_method"])
    instance = args["instance"]
    seed = args["seed"]
    threads = args["threads"]
    time_limit = args["time_limit"]
    gap_tolerance = args["gap_tolerance"]
    numeric_tolerance = args["numeric_tolerance"]
    result_path = output_path(
        args["output_dir"],
        instance,
        master_solver,
        sub_solver,
        sub_method,
        seed,
    )

    metadata = (
        timestamp_utc = Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sss") * "Z",
        problem = "uflp",
        dataset = "koerkel_ghosh",
        instance = instance,
        environment = "bnb",
        oracle = "classical",
        preprocessing = "none",
        master_solver = string(master_solver),
        sub_solver = string(sub_solver),
        sub_method = string(sub_method),
        threads = threads,
        julia_threads = Threads.nthreads(),
        seed = seed,
        time_limit = time_limit,
        gap_tolerance = gap_tolerance,
        numeric_tolerance = numeric_tolerance,
        warmup = !args["skip_warmup"],
        julia_version = string(VERSION),
        bendersx_version = string(Base.pkgversion(BendersX)),
    )

    @info "Starting one solver-composition run" instance master_solver sub_solver result_path
    try
        Random.seed!(seed)
        !args["skip_warmup"] && warmup!(
            master_solver,
            sub_solver,
            sub_method;
            threads = threads,
            seed = seed,
            gap_tolerance = gap_tolerance,
            numeric_tolerance = numeric_tolerance,
        )
        Random.seed!(seed)

        data = nothing
        env = nothing
        oracle = nothing
        build_time = @elapsed begin
            data = read_Simple_data(instance)
            env, oracle = build_environment(
                data,
                master_solver,
                sub_solver,
                sub_method;
                threads = threads,
                seed = seed,
                time_limit = time_limit,
                gap_tolerance = gap_tolerance,
                numeric_tolerance = numeric_tolerance,
            )
        end

        master_version = optimizer_version(env.master.model)
        sub_version = optimizer_version(oracle.model)
        summary = DataFrame()
        wall_time = @elapsed summary = solve!(env)
        result = summary[1, :]
        row = (
            metadata...,
            algorithm_status = string(nameof(typeof(env.termination_status))),
            build_time = build_time,
            wall_time = wall_time,
            reported_time = Float64(result.time),
            preprocessing_time = Float64(result.preprocessing_time),
            oracle_time = Float64(result.oracle_time),
            mean_oracle_time = Float64(result.mean_oracle_time),
            node_count = Int(result.node_count),
            callback_calls = Int(result.callback_calls),
            n_lazy_cuts = Int(result.n_lazy_cuts),
            n_user_cuts = Int(result.n_user_cuts),
            obj_bound = Float64(result.obj_bound),
            obj_val = Float64(result.obj_val),
            rel_gap = Float64(result.rel_gap),
            master_solver_version = master_version,
            sub_solver_version = sub_version,
            error = "",
        )
        write_result(result_path, row; overwrite = args["overwrite"])
        @info "Solver-composition run finished" status = row.algorithm_status obj_val = row.obj_val wall_time result_path
        return nothing
    catch error
        row = (
            metadata...,
            algorithm_status = "Error",
            build_time = missing,
            wall_time = missing,
            reported_time = missing,
            preprocessing_time = missing,
            oracle_time = missing,
            mean_oracle_time = missing,
            node_count = missing,
            callback_calls = missing,
            n_lazy_cuts = missing,
            n_user_cuts = missing,
            obj_bound = missing,
            obj_val = missing,
            rel_gap = missing,
            master_solver_version = "unknown",
            sub_solver_version = "unknown",
            error = sprint(showerror, error),
        )
        write_result(result_path, row; overwrite = args["overwrite"])
        rethrow()
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
