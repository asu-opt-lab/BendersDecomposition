export create_master_problem, create_sub_problem, create_dcglp
"""
    create_master_problem(data::AbstractData, cut_strategy::CutStrategy
)

Create the master problem formulation for Benders decomposition.

# Arguments
- `data::AbstractData`: Problem-specific data structure
- `cut_strategy::CutStrategy
`: Strategy for generating Benders cuts

# Returns
A concrete subtype of `AbstractMasterProblem` 

# Throws
- `MethodError`: If no implementation exists for the given data and cut strategy types
"""
function create_master_problem end

"""
    create_sub_problem(data::AbstractData, cut_strategy::CutStrategy
)

Create the sub problem formulation for Benders decomposition.

# Arguments
- `data::AbstractData`: Problem-specific data structure
- `cut_strategy::CutStrategy
`: Strategy for generating Benders cuts

# Returns
A concrete subtype of `AbstractSubProblem` 

# Throws
- `MethodError`: If no implementation exists for the given data and cut strategy types
"""
function create_sub_problem end

"""
    create_dcglp(data::AbstractData, cut_strategy::CutStrategy
)

Create the Dual Cut Generation Linear Program (DCGLP) formulation.

# Arguments
- `data::AbstractData`: Problem-specific data structure
- `cut_strategy::CutStrategy
`: Strategy for generating Benders cuts

# Returns
`DCGLP` structure 

# Throws
- `MethodError`: If no implementation exists for the given data and cut strategy types
"""
function create_dcglp end


function create_master_problem(data::AbstractData, cut_strategy::CutStrategy)
    throw(MethodError(create_master_problem, (data, cut_strategy)))
end

function create_sub_problem(data::AbstractData, cut_strategy::CutStrategy)
    throw(MethodError(create_sub_problem, (data, cut_strategy)))
end

function create_dcglp(data::AbstractData, cut_strategy::CutStrategy)
    throw(MethodError(create_dcglp, (data, cut_strategy)))
end

"""
Helper functions for DCGLP creation and manipulation
"""

function create_master_problem(data::AbstractData, cut_strategy::DisjunctiveCut)
    master = create_master_problem(data, cut_strategy.base_cut_strategy)
    relax_integrality(master.model)
    return master
end

function create_t_variable(model::Model, ::CutStrategy, data::AbstractData) 
    @variable(model, t >= -1e6)
end 


function create_sub_problem(data::AbstractData, cut_strategy::DisjunctiveCut)
    sub = create_sub_problem(data, cut_strategy.base_cut_strategy)
    return sub
end

# Include model-specific implementations
include("CFLP/cflp.jl")
include("UFLP/uflp.jl")
include("SCFLP/scflp.jl")
include("MCNDP/mcndp.jl")
include("SNIP/snip.jl")
include("base_dcglp.jl")