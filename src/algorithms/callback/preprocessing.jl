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
- `seq_type::Type{<:AbstractBendersSeq}`: Type of Benders sequence to use
- `params::AbstractBendersSeqParam`: Parameters for the Benders sequence
"""
struct RootNodePreprocessing <: AbstractRootNodePreprocessing
    oracle::AbstractOracle
    seq_type::Type{<:AbstractBendersSeq}
    params::AbstractBendersSeqParam

    function RootNodePreprocessing(oracle::AbstractOracle, seq_type::Type{<:AbstractBendersSeq}, params::AbstractBendersSeqParam)
        new(oracle, seq_type, params)
    end

    function RootNodePreprocessing(oracle::AbstractOracle; params::AbstractBendersSeqParam = BendersSeqParam())
        new(oracle, BendersSeq, params)
    end

    function RootNodePreprocessing(data::Data; params::AbstractBendersSeqParam = BendersSeqParam())
        new(ClassicalOracle(data), BendersSeq, params)
    end
end

"""
    root_node_processing!(data::Data, master::AbstractMaster, root_preprocessing::RootNodePreprocessing)

Process the root node of the branch-and-bound tree by temporarily relaxing integrality 
constraints and generating initial Benders cuts.

# Arguments
- `data::Data`: Problem data
- `master::AbstractMaster`: Master problem
- `root_preprocessing::RootNodePreprocessing`: Configuration for root node preprocessing

# Returns
- `Float64`: Time taken for root node processing
"""
function root_node_processing!(data::Data, master::AbstractMaster, root_preprocessing::RootNodePreprocessing)
    root_param = deepcopy(root_preprocessing.params)

    undo = relax_integrality(master.model)
    
    root_node_time = @elapsed begin
        BendersRootSeq = root_preprocessing.seq_type(data, master, root_preprocessing.oracle; param=root_param)
        solve!(BendersRootSeq)
    end
    
    undo()
    
    return root_node_time
end