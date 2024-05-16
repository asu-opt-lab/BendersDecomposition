function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--instance"
        help = "instance"
        default = "f100-c200-r5.0-p1"
        arg_type = AbstractString
        "--cut_strategy"
        help = "cut_strategy"
        default = ORDINARY_CUTSTRATEGY
        arg_type = Union{AbstractCutStrategy, Nothing}
        "--SplitCGLPNormType"
        help = "SplitCGLPNormType"
        default = nothing
        arg_type = Union{AbstractNormType, Nothing}
        "--SplitSetSelectionPolicy"
        help = "SplitSetSelectionPolicy"
        default = nothing
        arg_type = Union{AbstractSplitSetSelectionPolicy, Nothing}
        "--StrengthenCutStrategy"
        help = "StrengthenCutStrategy"
        default = nothing
        arg_type = Union{AbstractSplitStengtheningPolicy, Nothing}
        "--SplitBendersStrategy"
        help = "SplitBendersStrategy"
        default = nothing
        arg_type = Union{AbstractSplitBendersPolicy, Nothing}

    end

    return parse_args(s)
end