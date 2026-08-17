"""
    LPRelaxationPreprocessing <: AbstractBendersPreprocessing

Root-node preprocessing using a sequential Benders environment.

Before branch-and-bound begins, this strategy temporarily relaxes the integrality constraints of the master problem and runs the configured sequential Benders environment to generate initial cuts.

# Fields
- `oracle::AbstractOracle`: Oracle used to generate Benders cuts.
- `seq_type::Type{<:AbstractBendersSeq}`: Sequential Benders environment type used for preprocessing.
- `param::AbstractBendersSeqParam`: Parameters passed to the preprocessing environment.
"""
mutable struct LPRelaxationPreprocessing <: AbstractBendersPreprocessing
    oracle::AbstractOracle
    seq_env_type::Type{<:AbstractBendersSeq}
    param::AbstractBendersSeqParam

    function LPRelaxationPreprocessing(oracle::AbstractOracle; seq_env_type::Type{<:AbstractBendersSeq} = BendersSeq, param::AbstractBendersSeqParam = BendersSeqParam())
        new(oracle, seq_env_type, param)
    end
end

"""
    preprocess!(master::AbstractMaster, preprocessing::RootNodePreprocessing)

Process the root node of the branch-and-bound tree by temporarily relaxing integrality 
constraints and generating initial Benders cuts.

# Arguments
- `master::AbstractMaster`: Master problem
- `preprocessing::RootNodePreprocessing`: Configuration for root node preprocessing

# Returns
- `Float64`: Time taken for root node processing
"""
function preprocess!(master::AbstractMaster, preprocessing::LPRelaxationPreprocessing)
    # Use a private copy because preprocessing consumes part of the time budget.
    preprocessing_seq_param = deepcopy(preprocessing.param)

    # Relax integrality, ensure undo() always runs even on error
    undo = relax_integrality(master.model)
    
    # measure time and ensure undo() is called even if solve! errors
    tic = time()
    try
        env_preprocessing = preprocessing.seq_env_type(master, preprocessing.oracle; param = preprocessing_seq_param)
        solve!(env_preprocessing)
    finally
        # always restore integrality (even on exceptions)
        undo()
    end
    
    return time() - tic
end