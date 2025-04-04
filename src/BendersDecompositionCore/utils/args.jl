using ArgParse
using TOML
export parse_commandline, load_benders_params, load_cut_strategy, load_all_you_need, load_snip_data

function parse_commandline(;is_snip=false)
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--instance"
            help = "Instance name (overrides config file)"
            arg_type = is_snip ? Int : String
            required = true
            default = is_snip ? nothing : "f10-c10-r3-1"
        "--config"
            help = "Path to configuration file"
            default = "src/utils/config.toml"
            arg_type = String
        "--output_dir"
            help = "Output directory"
            default = "experiments"
            arg_type = String
        "--snip_no"
            help = "SNIP number"
            default = 1
            arg_type = Int
            required = false
        "--budget"
            help = "Budget"
            default = 30.0
            arg_type = Float64
            required = false
    end

    return parse_args(s)
end

function load_benders_params(config::Dict)
    alg_params = get(config, "algorithm_params", Dict())
    solver_params = get(config, "solver_params", Dict())
    cut_type = get(config["cut_strategy"], "type", "CLASSICAL_CUT")
    
    time_limit = get(alg_params, "time_limit", 3600)
    gap_tolerance = get(alg_params, "gap_tolerance", 1e-4)
    solver = get(alg_params, "solver", "GUROBI")
    verbose = get(alg_params, "verbose", true)
    
    return BendersParams(
        time_limit,
        gap_tolerance,
        solver,
        get(solver_params, "master", Dict()),
        get(solver_params, "sub", Dict()),
        cut_type == "DISJUNCTIVE_CUT" ? get(solver_params, "dcglp", Dict()) : Dict{String,Any}(),
        verbose
    )
end

const CUT_STRATEGY_MAP = Dict(
    "CLASSICAL_CUT" => ClassicalCut,
    "FAT_KNAPSACK_CUT" => FatKnapsackCut,
    "SLIM_KNAPSACK_CUT" => SlimKnapsackCut,
    "KNAPSACK_CUT" => KnapsackCut
)

const NORM_TYPE_MAP = Dict(
    "STANDARD_NORM" => StandardNorm,
    "L1NORM" => L1Norm,
    "L2NORM" => L2Norm,
    "LINFNORM" => LInfNorm
)

const STRENGTHENING_MAP = Dict(
    "PURE_DISJUNCTION" => PureDisjunctiveCut,
    "STRENGTHENED_DISJUNCTION" => StrengthenedDisjunctiveCut
)

function load_cut_strategy(config::Dict)
    cut_config = get(config, "cut_strategy", Dict())
    cut_type = get(cut_config, "type", "CLASSICAL_CUT")
    
    if cut_type == "DISJUNCTIVE_CUT"
        base_strategy = get(CUT_STRATEGY_MAP, cut_config["base_cut_strategy"]) do
            error("Unknown base cut strategy: $(cut_config["base_cut_strategy"])")
        end
        
        norm_type = get(NORM_TYPE_MAP, cut_config["norm_type"]) do
            error("Unknown norm type: $(cut_config["norm_type"])")
        end
        
        strengthening = get(STRENGTHENING_MAP, cut_config["cut_strengthening"]) do
            error("Unknown cut strengthening type: $(cut_config["cut_strengthening"])")
        end
        
        return DisjunctiveCut(
            base_strategy(),
            norm_type(),
            strengthening(),
            get(cut_config, "use_two_sided_cuts", false),
            get(cut_config, "include_master_cuts", true),
            get(cut_config, "reuse_dcglp", true),
            get(cut_config, "verbose", false)
        )
    else
        return get(CUT_STRATEGY_MAP, cut_type) do
            error("Unknown cut strategy type: $cut_type")
        end()
    end
end

function load_base_config(;is_snip=false)
    args = parse_commandline(is_snip=is_snip)
    config = TOML.parsefile(args["config"])
    return args, config
end

function load_all_you_need()
    args, config = load_base_config()
    return args["instance"], 
           args["output_dir"], 
           load_cut_strategy(config), 
           load_benders_params(config)
end

function load_snip_data()
    args, config = load_base_config(is_snip=true)
    return args["instance"], 
           args["snip_no"], 
           args["budget"], 
           args["output_dir"], 
           load_cut_strategy(config), 
           load_benders_params(config)
end
