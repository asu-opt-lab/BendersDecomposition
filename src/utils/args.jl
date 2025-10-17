using ArgParse
export parse_commandline

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--instance"
            help = "Instance name (overrides config file)"
            default = "f10-c10-r3-1"
            arg_type = String
            required = true
        "--output_dir"
            help = "Output directory"
            default = "experiments"
            arg_type = String
        "--seed"
            help = "Random seed"
            default = 1234
            arg_type = Int
    end

    return parse_args(s)
end

