include("common.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--output_dir"
            arg_type = String
            default = RAW_RESULTS_DIR
        "--time_limit"
            arg_type = Float64
            default = 600.0
        "--gap_tolerance"
            help = "BendersSeq gap tolerance in percent."
            arg_type = Float64
            default = 1e-6
        "--solver_threads"
            arg_type = Int
            default = 1
        "--problem"
            help = "Optional filter: uflp_correctness, cflp_correctness, or all."
            arg_type = String
            default = "all"
        "--instance"
            help = "Optional single-instance filter, e.g. p1. Use all for the full list."
            arg_type = String
            default = "all"
        "--verbose"
            action = :store_true
    end
    return ArgParse.parse_args(s)
end

function main()
    args = parse_commandline()
    files = default_output_files("01_correctness", args["output_dir"])
    cases = [
        ("uflp_correctness", UFLP_CORRECTNESS_INSTANCES, ["classical", "unified", "pareto", "ufl_knapsack"]),
        ("cflp_correctness", CFLP_CORRECTNESS_INSTANCES, ["classical", "unified", "pareto", "cfl_knapsack"]),
    ]

    for (problem, instances, oracles) in cases
        args["problem"] == "all" || args["problem"] == problem || continue
        selected_instances = args["instance"] == "all" ? instances : [args["instance"]]
        for instance in instances
            instance in selected_instances || continue
            data, data_source = load_paper_instance(problem, instance)
            @info "baseline" problem instance
            baseline = log_baseline_case(;
                experiment = "01_correctness",
                problem = problem,
                instance = instance,
                data = data,
                data_source = data_source,
                repeat = 1,
                time_limit = args["time_limit"],
                solver_threads = args["solver_threads"],
                summary_file = files.summary,
                verbose = args["verbose"],
            )

            for oracle_name in oracles
                @info "benders" problem instance oracle_name
                run_benders_case(;
                    experiment = "01_correctness",
                    problem = problem,
                    instance = instance,
                    data = data,
                    data_source = data_source,
                    env_name = "seq",
                    oracle_name = oracle_name,
                    repeat = 1,
                    time_limit = args["time_limit"],
                    gap_tolerance = args["gap_tolerance"],
                    solver_threads = args["solver_threads"],
                    summary_file = files.summary,
                    trace_file = files.trace,
                    baseline_obj = baseline.obj_value,
                    verbose = args["verbose"],
                )
            end
        end
    end
end

main()
