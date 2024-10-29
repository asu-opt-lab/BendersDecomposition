using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--problem"
            help = "problem"
            default = "UFLP"
            arg_type = String
        "--instance"
            help = "instance"
            default = "f100-c100-r5.0-p1"
            arg_type = String
        "--output_dir"
            help = "output directory"
            default = "results"
            arg_type = String
        "--solver"
            help = "Solver selection"
            default = "Gurobi"
            arg_type = String
        "--cut_strategy"
            help = "Cut strategy"
            default = "STANDARD_CUTSTRATEGY"
            arg_type = String
    end

    parsed_args = parse_args(s)

    if parsed_args["cut_strategy"] == "SPLIT_CUTSTRATEGY"
        @add_arg_table! s begin
            "--norm_type"
                help = "Norm type"
                default = "L2GAMMANORM"
                arg_type = String
            "--split_set_selection"
                help = "Split set selection policy"
                default = "MOST_FRAC_INDEX"
                arg_type = String
            "--split_strengthening"
                help = "Split strengthening policy"
                default = "SPLIT_PURE_CUT_STRATEGY"
                arg_type = String
            "--split_benders"
                help = "Split Benders policy"
                default = "ALL_SPLIT_BENDERS_STRATEGY"
                arg_type = String
        end
        parsed_args = parse_args(s)
    end

    return parsed_args
end

args = parse_commandline()
