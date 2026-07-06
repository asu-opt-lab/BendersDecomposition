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
            help = "BendersSeq gap tolerance in percent."
            arg_type = Float64
            default = 1e-4
        "--solver_threads"
            arg_type = Int
            default = 1
        "--repeats"
            arg_type = Int
            default = 3
        "--repeat_index"
            help = "Optional single repeat. Use 0 to run 1:repeats."
            arg_type = Int
            default = 0
        "--problem"
            help = "Optional filter: uflp, cflp, scflp, or all."
            arg_type = String
            default = "all"
        "--instance"
            help = "Optional single-instance filter. Use all for the full list."
            arg_type = String
            default = "all"
        "--oracle"
            help = "Optional oracle filter. Use all for the experiment matrix."
            arg_type = String
            default = "all"
        "--verbose"
            action = :store_true
    end
    return ArgParse.parse_args(s)
end

function main()
    args = parse_commandline()
    files = default_output_files("02_oracle_ablation", args["output_dir"])
    cases = [
        ("uflp", UFLP_250_INSTANCES, ["classical", "unified", "pareto", "ufl_knapsack"]),
        ("cflp", CFLP_200_INSTANCES, ["classical", "unified", "pareto", "cfl_knapsack"]),
        ("scflp", SCFLP_100X200_INSTANCES, ["classical", "unified", "pareto", "cfl_knapsack"]),
    ]

    repeat_range = args["repeat_index"] == 0 ? (1:args["repeats"]) : (args["repeat_index"]:args["repeat_index"])
    for repeat in repeat_range
        for (problem, instances, oracles) in cases
            args["problem"] == "all" || args["problem"] == problem || continue
            selected_instances = args["instance"] == "all" ? instances : [args["instance"]]
            selected_oracles = args["oracle"] == "all" ? oracles : [args["oracle"]]
            for instance in instances
                instance in selected_instances || continue
                data, data_source = load_paper_instance(problem, instance)
                for oracle_name in selected_oracles
                    @info "oracle ablation" repeat problem instance oracle_name
                    run_benders_case(;
                        experiment = "02_oracle_ablation",
                        problem = problem,
                        instance = instance,
                        data = data,
                        data_source = data_source,
                        env_name = "seq",
                        oracle_name = oracle_name,
                        repeat = repeat,
                        time_limit = args["time_limit"],
                        gap_tolerance = args["gap_tolerance"],
                        solver_threads = args["solver_threads"],
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
