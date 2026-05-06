using ArgParse
export parse_commandline

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--instance"
            help = "Instance name (overrides config file)"
            default = "f10-c10-r3-1"
            arg_type = String
            required = false
        "--output_dir"
            help = "Output directory"
            default = "experiments"
            arg_type = String
        "--seed"
            help = "Random seed"
            default = 1234
            arg_type = Int
        "--snip_instance"
            help = "SNIP instance number"
            default = 0
            arg_type = Int
            required = false
        "--snip_no"
            help = "SNIP number"
            default = 1
            arg_type = Int
            required = false
        "--snip_budget"
            help = "SNIP budget"
            default = 30.0
            arg_type = Float64
            required = false
    end

    return parse_args(s)
end

