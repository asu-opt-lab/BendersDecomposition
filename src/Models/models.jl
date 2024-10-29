# include("CFLP/cflp.jl")
include("UFLP/uflp.jl")

export create_master_problem, create_sub_problem, create_dcglp


function create_master_problem(data::AbstractData, cut_strategy::CutGenerationStrategy)
    error("Unsupported cut strategy: $(typeof(cut_strategy))")
end

function create_sub_problem(data::AbstractData, cut_strategy::CutGenerationStrategy)
    error("Unsupported cut strategy: $(typeof(cut_strategy))")
end

function create_dcglp(data::AbstractData, disjunctive_inequality::Tuple{Vector{Float64}, Float64}, _cut_strategy::CutGenerationStrategy, norm_type::AbstractNormType)
    error("Unsupported cut strategy: $(typeof(_cut_strategy))")
end

