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
            help = "Environment for the SCFLP scaling run: callback, bnb, seq, or seq_inout."
            arg_type = String
            default = "callback"
        "--verbose"
            action = :store_true
    end
    return ArgParse.parse_args(s)
end

function main()
    args = parse_commandline()
    files = default_output_files("04_parallel_scflp", args["output_dir"])
    problem = "scflp"

    repeat_range = args["repeat_index"] == 0 ? (1:args["repeats"]) : (args["repeat_index"]:args["repeat_index"])
    selected_instances = args["instance"] == "all" ? SCFLP_50X100_INSTANCES : [args["instance"]]
    for repeat in repeat_range
        for instance in selected_instances
            data, data_source = load_paper_instance(problem, instance)
            @info "parallel scflp" repeat instance Threads.nthreads() args["oracle"] args["env"]
            run_benders_case(;
                experiment = "04_parallel_scflp",
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
            )
        end
    end
end

main()
