using ArgParse

function parse_commandline(problem_type::SNIP)
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--snip_instance"
            help = "SNIP instance number"
            default = 0
            arg_type = Int
            required = true
        "--snip_no"
            help = "SNIP number"
            default = 1
            arg_type = Int
            required = true
        "--snip_budget"
            help = "SNIP budget"
            default = 30.0
            arg_type = Float64
            required = true
        "--output_dir"
            help = "Output directory"
            default = "experiments"
            arg_type = String
    end

    return parse_args(s)
end

