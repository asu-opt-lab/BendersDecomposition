using ArgParse
using TOML
export parse_commandline, load_benders_params, load_cut_strategy, load_all_you_need

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--instance"
            help = "Instance name (overrides config file)"
            default = "f10-c10-r3-1"
            arg_type = String
        "--config"
            help = "Path to configuration file"
            default = "src/utils/config.toml"
            arg_type = String
        "--output_dir"
            help = "Output directory"
            default = "experiments"
            arg_type = String
    end

    return parse_args(s)
end


function load_benders_params(config)
    alg_params = config["algorithm_params"]
    solver_params = config["solver_params"]
    # Check if using disjunctive cuts
    is_disjunctive = get(config["cut_strategy"], "type", "STANDARD_CUT") == "DISJUNCTIVE_CUT"
    
    # Create BendersParams with values from config
    if is_disjunctive
        return BendersParams(
            alg_params["time_limit"],
            alg_params["gap_tolerance"],
            alg_params["solver"],
            solver_params["master"],
            solver_params["sub"],
            solver_params["dcglp"],  # Only for disjunctive cuts
            alg_params["verbose"],
        )
    else
        return BendersParams(
            alg_params["time_limit"],
            alg_params["gap_tolerance"],
            alg_params["solver"],
            solver_params["master"],
            solver_params["sub"],
            Dict{String,Any}(),  # dcglp solver params
            alg_params["verbose"],
        )
    end
end

function load_cut_strategy(config)
    cut_config = config["cut_strategy"]
    
    # Check cut strategy type
    if cut_config["type"] == "DISJUNCTIVE_CUT"
        # Convert string to corresponding cut strategy type
        base_strategy = if cut_config["base_cut_strategy"] == "STANDARD_CUT"
            ClassicalCut()
        elseif cut_config["base_cut_strategy"] == "FAT_KNAPSACK_CUT"
            FatKnapsackCut()
        elseif cut_config["base_cut_strategy"] == "SLIM_KNAPSACK_CUT"
            SlimKnapsackCut()
        elseif cut_config["base_cut_strategy"] == "KNAPSACK_CUT"
            KnapsackCut()
        else
            error("Unknown base cut strategy: $(cut_config["base_cut_strategy"])")
        end
        
        # Convert string to corresponding norm type
        norm_type = if cut_config["norm_type"] == "STANDARD_NORM"
            StandardNorm()
        elseif cut_config["norm_type"] == "L1NORM"
            L1Norm()
        elseif cut_config["norm_type"] == "L2NORM"
            L2Norm()
        elseif cut_config["norm_type"] == "LINFNORM"
            LInfNorm()
        else
            error("Unknown norm type: $(cut_config["norm_type"])")
        end
        
        # Convert string to corresponding cut strengthening type
        strengthening = if cut_config["cut_strengthening"] == "PURE_DISJUNCTION"
            PureDisjunctiveCut()
        elseif cut_config["cut_strengthening"] == "STRENGTHENED_DISJUNCTION"
            StrengthenedDisjunctiveCut()
        else
            error("Unknown cut strengthening type: $(cut_config["cut_strengthening"])")
        end
        
        # Create and return the DisjunctiveCut


        return DisjunctiveCut(
            base_strategy,
            norm_type,
            strengthening,
            cut_config["use_two_sided_cuts"],
            cut_config["include_master_cuts"],
            cut_config["reuse_dcglp"],
            cut_config["verbose"]
            )
    else
        # Handle standard cut strategy
        return if cut_config["type"] == "STANDARD_CUT"
            ClassicalCut()
        elseif cut_config["type"] == "KNAPSACK_CUT"
            KnapsackCut()
        else
            error("Unknown cut strategy type: $(cut_config["type"])")
        end
    end
end

# Update the main configuration loading
function load_all_you_need()
    args = parse_commandline()
    instance = args["instance"]
    output_dir = args["output_dir"]
    config = TOML.parsefile(args["config"])
    @info config
    
    benders_params = load_benders_params(config)
    cut_strategy = load_cut_strategy(config)
    
    return instance, output_dir, cut_strategy, benders_params
end



