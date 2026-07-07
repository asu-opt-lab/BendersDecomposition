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
            help = "Keep this at 1 for parallel separation experiments."
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
        "--instance"
            help = "Optional single-instance filter. Use all for the full list."
            arg_type = String
            default = "all"
        "--oracle"
            arg_type = String
            default = "cfl_knapsack"
        "--env"
            help = "Environment for the LP-relaxed SCFLP run: seq or seq_inout."
            arg_type = String
            default = "seq"
        "--verbose"
            action = :store_true
    end
    return ArgParse.parse_args(s)
end

function main()
    args = parse_commandline()
    files = default_output_files("05_parallel_scflp_lp_relax", args["output_dir"])
    problem = "scflp"

    if !(args["env"] in ("seq", "seq_inout"))
        error("05_parallel_scflp_lp_relax expects --env seq or --env seq_inout, got $(args["env"])")
    end

    repeat_range = args["repeat_index"] == 0 ? (1:args["repeats"]) : (args["repeat_index"]:args["repeat_index"])
    selected_instances = args["instance"] == "all" ? SCFLP_50X50_S512_INSTANCES : [args["instance"]]
    for repeat in repeat_range
        for instance in selected_instances
            data, data_source = load_paper_instance(problem, instance)
            @info "parallel scflp lp relaxation" repeat instance Threads.nthreads() args["oracle"] args["env"]
            run_benders_case(;
                experiment = "05_parallel_scflp_lp_relax",
                problem = problem,
                instance = instance,
                data = data,
                data_source = data_source,
                env_name = args["env"],
                oracle_name = args["oracle"],
                repeat = repeat,
                time_limit = args["time_limit"],
                gap_tolerance = args["gap_tolerance"],
                solver_threads = args["solver_threads"],
                solver_name = args["solver"],
                summary_file = files.summary,
                trace_file = files.trace,
                verbose = args["verbose"],
                lp_relaxation = true,
            )
        end
    end
end

main()
