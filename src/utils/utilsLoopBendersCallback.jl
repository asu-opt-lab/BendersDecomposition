export BendersCallbackParam

abstract type AbstractBendersCallbackParam <: AbstractLoopParam end

"""
Parameters for configuring the Callback-based Benders decomposition algorithm.

Contains settings for:
- `time_limit`: Maximum runtime allowed for the algorithm in seconds.
- `gap_tolerance`: Relative optimality gap tolerance for termination.
- `verbose`: Controls the level of logging output during execution.
- `preprocessing_type`: Type of preprocessing to apply at the root node (NoPreprocessing, SeqPreprocessing, or SeqInOutPreprocessing).
- `root_param`: Parameters for root node preprocessing.

These parameters allow fine-tuning of the Benders algorithm performance.
"""
# mutable struct BendersCallbackParam <: AbstractBendersCallbackParam
#     time_limit::Float64
#     gap_tolerance::Float64
#     verbose::Bool
#     preprocessing_type::RootNodePreprocessingType
#     root_param::Union{Nothing,AbstractBendersSeqParam}

#     function BendersCallbackParam(; 
#                         time_limit::Float64 = 7200.0, 
#                         gap_tolerance::Float64 = 1e-6, 
#                         verbose::Bool = true,
#                         preprocessing_type::RootNodePreprocessingType = NoPreprocessing(),
#                         root_param::Union{Nothing,AbstractBendersSeqParam} = nothing
#                         ) 
        
#         new(time_limit, gap_tolerance, verbose, 
#             preprocessing_type, root_param)
#     end
# end 

mutable struct BendersCallbackParam <: AbstractBendersCallbackParam
    time_limit::Float64
    gap_tolerance::Float64
    verbose::Bool
    preprocessing_type::Any
    root_param::Union{Nothing,AbstractBendersSeqParam}

    function BendersCallbackParam(; 
                        time_limit::Float64 = 7200.0, 
                        gap_tolerance::Float64 = 1e-6, 
                        verbose::Bool = true,
                        preprocessing_type::Any = nothing,
                        root_param::Union{Nothing,AbstractBendersSeqParam} = nothing
                        ) 
        
        new(time_limit, gap_tolerance, verbose, 
            preprocessing_type, root_param)
    end
end 