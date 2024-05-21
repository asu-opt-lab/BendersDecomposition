function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--instance"
        help = "instance"
        default = "f100-c200-r5.0-p1"
        arg_type = AbstractString
        "--cut_strategy"
        help = "cut_strategy"
        default = "ORDINARY_CUTSTRATEGY"
        arg_type = AbstractString
        "--SplitCGLPNormType"
        help = "SplitCGLPNormType"
        default = "nothing"
        arg_type = AbstractString
        "--SplitSetSelectionPolicy"
        help = "SplitSetSelectionPolicy"
        default = "nothing"
        arg_type = AbstractString
        "--StrengthenCutStrategy"
        help = "StrengthenCutStrategy"
        default = "nothing"
        arg_type = AbstractString
        "--SplitBendersStrategy"
        help = "SplitBendersStrategy"
        default = "nothing"
        arg_type = AbstractString

    end

    return parse_args(s)
end

function set_params_attribute(algo_params,::Type{AbstractCutStrategy},s::AbstractString)
    if s == "ORDINARY_CUTSTRATEGY" 
        algo_params.cut_strategy = ORDINARY_CUTSTRATEGY
    elseif s == "SPLIT_CUTSTRATEGY"
        algo_params.cut_strategy = SPLIT_CUTSTRATEGY
    elseif s == "nothing"
        algo_params.cut_strategy = nothing
    else
        throw(ArgumentError("Invalid cut strategy"))
    end
end

function set_params_attribute(algo_params,::Type{AbstractNormType},s::AbstractString)
    if s == "STANDARDNORM"
        algo_params.SplitCGLPNormType = STANDARDNORM
    elseif s == "L1GammaNorm"
        algo_params.SplitCGLPNormType = L1GAMMANORM
    elseif s == "L2GammaNorm"
        algo_params.SplitCGLPNormType = L2GAMMANORM
    elseif s == "LINFGAMMANORM"
        algo_params.SplitCGLPNormType = LINFGAMMANORM
    elseif s == "nothing"
        algo_params.SplitCGLPNormType = nothing
    else
        throw(ArgumentError("Invalid norm type"))
    end
end

function set_params_attribute(algo_params,::Type{AbstractSplitSetSelectionPolicy},s::AbstractString)
    if s == "MOST_FRAC_INDEX"
        algo_params.SplitSetSelectionPolicy = MOST_FRAC_INDEX
    elseif s == "RANDOM_INDEX"
        algo_params.SplitSetSelectionPolicy = RANDOM_INDEX
    elseif s == "nothing"
        algo_params.SplitSetSelectionPolicy = nothing
    else
        throw(ArgumentError("Invalid split set selection policy"))
    end
end

function set_params_attribute(algo_params,::Type{AbstractSplitBendersPolicy},s::AbstractString)
    if s == "NO_SPLIT_BENDERS_STRATEGY"
        algo_params.SplitBendersStrategy = NO_SPLIT_BENDERS_STRATEGY
    elseif s == "ALL_SPLIT_BENDERS_STRATEGY"
        algo_params.SplitBendersStrategy = ALL_SPLIT_BENDERS_STRATEGY
    elseif s == "TIGHT_SPLIT_BENDERS_STRATEGY"
        algo_params.SplitBendersStrategy = TIGHT_SPLIT_BENDERS_STRATEGY
    elseif s == "nothing"
        algo_params.SplitBendersStrategy = nothing
    else
        throw(ArgumentError("Invalid split benders policy"))
    end
end

function set_params_attribute(algo_params,::Type{AbstractSplitStengtheningPolicy},s::AbstractString)
    if s == "SPLIT_PURE_CUT_STRATEGY"
        algo_params.StrengthenCutStrategy = SPLIT_PURE_CUT_STRATEGY
    elseif s == "SPLIT_STRENGTHEN_CUT_STRATEGY"
        algo_params.StrengthenCutStrategy = SPLIT_STRENGTHEN_CUT_STRATEGY
    elseif s == "nothing"
        algo_params.StrengthenCutStrategy = nothing
    else
        throw(ArgumentError("Invalid split strengthening policy"))
    end
end































# function parse_commandline()
#     s = ArgParseSettings()

#     @add_arg_table s begin
#         "--instance"
#         help = "instance"
#         default = "f100-c200-r5.0-p1"
#         arg_type = AbstractString
#         "--cut_strategy"
#         help = "cut_strategy"
#         default = ORDINARY_CUTSTRATEGY
#         arg_type = Union{AbstractCutStrategy, Nothing}
#         "--SplitCGLPNormType"
#         help = "SplitCGLPNormType"
#         default = nothing
#         arg_type = Union{AbstractNormType, Nothing}
#         "--SplitSetSelectionPolicy"
#         help = "SplitSetSelectionPolicy"
#         default = nothing
#         arg_type = Union{AbstractSplitSetSelectionPolicy, Nothing}
#         "--StrengthenCutStrategy"
#         help = "StrengthenCutStrategy"
#         default = nothing
#         arg_type = Union{AbstractSplitStengtheningPolicy, Nothing}
#         "--SplitBendersStrategy"
#         help = "SplitBendersStrategy"
#         default = nothing
#         arg_type = Union{AbstractSplitBendersPolicy, Nothing}

#     end

#     return parse_args(s)
# end

# function ArgParse.parse_item(::Type{Union{AbstractCutStrategy, Nothing}}, s::AbstractString)
#     if s == "ORDINARY_CUTSTRATEGY" 
#         return ORDINARY_CUTSTRATEGY
#     elseif s == "SPLIT_CUTSTRATEGY"
#         return SPLIT_CUTSTRATEGY
#     elseif isnothing(s)
#         return "nothing"
#     else
#         @info s
#         throw(ArgumentError("Invalid cut strategy"))
#     end
# end

# function reConvert_cutstrategy(s::Union{AbstractCutStrategy, Nothing})
#     if s == ORDINARY_CUTSTRATEGY
#         return "ORDINARY_CUTSTRATEGY"
#     elseif s == SPLIT_CUTSTRATEGY
#         return "SPLIT_CUTSTRATEGY"
#     elseif isnothing(s)
#         return "nothing"
#     else
#         throw(ArgumentError("Invalid cut strategy"))
#     end
# end


# function ArgParse.parse_item(::Type{Union{AbstractNormType, Nothing}}, s::AbstractString)
#     if s == "STANDARDNORM"
#         return STANDARDNORM
#     elseif s == "L1GammaNorm"
#         return L1GammaNorm
#     elseif s == "L2GammaNorm"
#         return L2GammaNorm
#     elseif s == "LINFGAMMANORM"
#         return LINFGAMMANORM
#     else
#         throw(ArgumentError("Invalid norm type"))
#     end
# end

# function reConvert_normtype(s::Union{AbstractNormType, Nothing})
#     if s == STANDARDNORM
#         return "STANDARDNORM"
#     elseif s == L1GAMMANORM
#         return "L1GammaNorm"
#     elseif s == L2GAMMANORM
#         return "L2GammaNorm"
#     elseif s == LINFGAMMANORM
#         return "LINFGAMMANORM"
#     elseif isnothing(s)
#         return "nothing"
#     else
#         throw(ArgumentError("Invalid norm type"))
#     end
# end

# function ArgParse.parse_item(::Type{Union{AbstractSplitSetSelectionPolicy, Nothing}}, s::AbstractString)
#     if s == "MOST_FRAC_INDEX"
#         return MOST_FRAC_INDEX
#     elseif s == "RANDOM_INDEX"
#         return RANDOM_INDEX
#     else
#         throw(ArgumentError("Invalid split set selection policy"))
#     end
# end

# function reConvert_splitsetselectionpolicy(s::Union{AbstractSplitSetSelectionPolicy, Nothing})
#     if s == MOST_FRAC_INDEX
#         return "MOST_FRAC_INDEX"
#     elseif s == RANDOM_INDEX
#         return "RANDOM_INDEX"
#     elseif isnothing(s)
#         return "nothing"
#     else
#         throw(ArgumentError("Invalid split set selection policy"))
#     end
# end

# function ArgParse.parse_item(::Type{Union{AbstractSplitBendersPolicy, Nothing}}, s::AbstractString)
#     if s == "NO_SPLIT_BENDERS_STRATEGY"
#         return NO_SPLIT_BENDERS_STRATEGY
#     elseif s == "ALL_SPLIT_BENDERS_STRATEGY"
#         return ALL_SPLIT_BENDERS_STRATEGY
#     elseif s == "TIGHT_SPLIT_BENDERS_STRATEGY"
#         return TIGHT_SPLIT_BENDERS_STRATEGY
#     else
#         throw(ArgumentError("Invalid split benders policy"))
#     end
# end

# function reConvert_splitbenderspolicy(s::Union{AbstractSplitBendersPolicy, Nothing})
#     if s == NO_SPLIT_BENDERS_STRATEGY
#         return "NO_SPLIT_BENDERS_STRATEGY"
#     elseif s == ALL_SPLIT_BENDERS_STRATEGY
#         return "ALL_SPLIT_BENDERS_STRATEGY"
#     elseif s == TIGHT_SPLIT_BENDERS_STRATEGY
#         return "TIGHT_SPLIT_BENDERS_STRATEGY"
#     elseif isnothing(s)
#         return "nothing"
#     else
#         throw(ArgumentError("Invalid split benders policy"))
#     end
# end

# function ArgParse.parse_item(::Type{Union{AbstractSplitStengtheningPolicy, Nothing}}, s::AbstractString)
#     if s == "SPLIT_PURE_CUT_STRATEGY"
#         return SPLIT_PURE_CUT_STRATEGY
#     elseif s == "SPLIT_STRENGTHEN_CUT_STRATEGY"
#         return SPLIT_STRENGTHEN_CUT_STRATEGY
#     else
#         throw(ArgumentError("Invalid split strengthening policy"))
#     end
# end

# function reConvert_splitstengtheningpolicy(s::Union{AbstractSplitStengtheningPolicy, Nothing})
#     if s == SPLIT_PURE_CUT_STRATEGY
#         return "SPLIT_PURE_CUT_STRATEGY"
#     elseif s == SPLIT_STRENGTHEN_CUT_STRATEGY
#         return "SPLIT_STRENGTHEN_CUT_STRATEGY"
#     elseif isnothing(s)
#         return "nothing"
#     else
#         throw(ArgumentError("Invalid split strengthening policy"))
#     end
# end