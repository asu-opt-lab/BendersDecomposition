include("common.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--output_dir"
            arg_type = String
            default = RAW_RESULTS_DIR
        "--time_limit"
            arg_type = Float64
            default = 1800.0
        "--gap_tolerance"
            help = "Gap tolerance passed to the selected environment."
            arg_type = Float64
            default = 1e-4
        "--solver_threads"
            arg_type = Int
            default = 1
        "--solver"
            help = "MILP/LP solver backend: gurobi or cplex."
            arg_type = String
            default = DEFAULT_SOLVER
        "--repeats"
            arg_type = Int
            default = 3
        "--repeat_index"
            help = "Optional single repeat. Use 0 to run 1:repeats."
            arg_type = Int
            default = 0
        "--problem"
            help = "Optional filter: uflp or all."
            arg_type = String
            default = "all"
        "--instance"
            help = "Optional single-instance filter. Use all for the full list."
            arg_type = String
            default = "all"
        "--env"
            help = "Optional environment filter: seq, seq_inout, callback, or all."
            arg_type = String
            default = "all"
        "--verbose"
            action = :store_true
    end
    return ArgParse.parse_args(s)
end

function main()
    args = parse_commandline()
    files = default_output_files("03_environment_ablation", args["output_dir"])
    cases = [
        ("uflp", UFLP_250_INSTANCES, "ufl_knapsack"),
    ]
    envs = ["seq", "seq_inout", "callback"]

    repeat_range = args["repeat_index"] == 0 ? (1:args["repeats"]) : (args["repeat_index"]:args["repeat_index"])
    for repeat in repeat_range
        for (problem, instances, oracle_name) in cases
            args["problem"] == "all" || args["problem"] == problem || continue
            if args["instance"] != "all" && !(args["instance"] in instances)
                error("Unknown $(problem) instance $(args["instance"]). Expected one of: $(join(instances, ", ")).")
            end
            if args["env"] != "all" && !(args["env"] in envs)
                error("Unknown environment $(args["env"]). Expected one of: $(join(envs, ", ")).")
            end
            selected_instances = args["instance"] == "all" ? instances : [args["instance"]]
            selected_envs = args["env"] == "all" ? envs : [args["env"]]
            for instance in instances
                instance in selected_instances || continue
                data, data_source = load_paper_instance(problem, instance)
                for env_name in selected_envs
                    @info "environment ablation" repeat problem instance env_name oracle_name
                    run_benders_case(;
                        experiment = "03_environment_ablation",
                        problem = problem,
                        instance = instance,
                        data = data,
                        data_source = data_source,
                        env_name = env_name,
                        oracle_name = oracle_name,
                        repeat = repeat,
                        time_limit = args["time_limit"],
                        gap_tolerance = args["gap_tolerance"],
                        solver_threads = args["solver_threads"],
                        solver_name = args["solver"],
                        summary_file = files.summary,
                        trace_file = files.trace,
                        verbose = args["verbose"],
                    )
                end
            end
        end
    end
end

main()
