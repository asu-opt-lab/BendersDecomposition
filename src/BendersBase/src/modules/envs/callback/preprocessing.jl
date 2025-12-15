export RootNodePreprocessing, NoRootNodePreprocessing
"""
    AbstractRootNodePreprocessing

Abstract type for root node preprocessing in Benders decomposition.
"""
abstract type AbstractRootNodePreprocessing end

"""
    NoRootNodePreprocessing <: AbstractRootNodePreprocessing

Indicates that no preprocessing should be done at the root node of the branch-and-bound tree.
"""
struct NoRootNodePreprocessing <: AbstractRootNodePreprocessing end

"""
    RootNodePreprocessing <: AbstractRootNodePreprocessing

Represents preprocessing to be performed at the root node of the branch-and-bound tree.
Used to generate initial cuts before the branch-and-bound procedure begins.

# Fields
- `oracle::AbstractOracle`: Oracle used to generate Benders cuts
- `seq_type::Type{<:AbstractBendersSeq}`: Type of BendersSeq to use
- `params::AbstractBendersSeqParam`: Parameters for the BendersSeq
"""
mutable struct RootNodePreprocessing <: AbstractRootNodePreprocessing
    oracle::AbstractOracle
    seq_type::Type{<:AbstractBendersSeq}
    params::AbstractBendersSeqParam

    function RootNodePreprocessing(oracle::AbstractOracle, seq_type::Type{<:AbstractBendersSeq}, params::AbstractBendersSeqParam)
        new(oracle, seq_type, params)
    end

    function RootNodePreprocessing(oracle::AbstractOracle; params::AbstractBendersSeqParam = BendersSeqParam())
        new(oracle, BendersSeq, params)
    end
end

"""
    root_node_processing!(master::AbstractMaster, root_preprocessing::NoRootNodePreprocessing) -> Float64
    
No-op implementation for NoRootNodePreprocessing.
# Returns
- `Float64`: 0.0 (no time spent)
"""
function root_node_processing!(master::AbstractMaster, root_preprocessing::NoRootNodePreprocessing)
    return 0.0
end

"""
    root_node_processing!(master::AbstractMaster, root_preprocessing::RootNodePreprocessing)

Process the root node of the branch-and-bound tree by temporarily relaxing integrality 
constraints and generating initial Benders cuts.

# Arguments
- `master::AbstractMaster`: Master problem
- `root_preprocessing::RootNodePreprocessing`: Configuration for root node preprocessing

# Returns
- `Float64`: Time taken for root node processing
"""
function root_node_processing!(master::AbstractMaster, root_preprocessing::RootNodePreprocessing)
    root_param = deepcopy(root_preprocessing.params)

    # Relax integrality, ensure undo() always runs even on error
    undo = relax_integrality(master.model)
    
    # measure time and ensure undo() is called even if solve! errors
    tic = time()
    try
        env_root = root_preprocessing.seq_type(master, root_preprocessing.oracle; param = root_param)
        solve!(env_root)
    finally
        # always restore integrality (even on exceptions)
        undo()
    end
    
    return time() - tic
end