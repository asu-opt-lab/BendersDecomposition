export DisjunctiveRootNodePreprocessing

"""
    DisjunctiveRootNodePreprocessing <: AbstractRootNodePreprocessing

Two-phase root-node preprocessing that first runs a typical
oracle-based separation and then a disjunctive oracle separation to
produce stronger initial cuts.
# Fields
- `typical_oracle::AbstractTypicalOracle`: Oracle for the first phase (classical preprocessing)
- `disjunctive_oracle::AbstractDisjunctiveOracle`: Oracle for the second phase
- `seq_type::Type{<:AbstractBendersSeq}`: Type of sequential Benders algorithm to use
- `params::AbstractBendersSeqParam`: Parameters for the sequential algorithm

# Constructor
```julia
DisjunctiveRootNodePreprocessing(
    typical_oracle::AbstractTypicalOracle,
    disjunctive_oracle::AbstractDisjunctiveOracle;
    seq_type::Type{<:AbstractBendersSeq} = BendersSeq,
    params::AbstractBendersSeqParam = BendersSeqParam()
)
```

# Examples
```julia
preprocessing = DisjunctiveRootNodePreprocessing(typical_oracle, disj_oracle)
# Use with BendersBnB
env = BendersBnB(master, preprocessing, lazy_callback, user_callback)
```
See also: [`RootNodePreprocessing`](@ref), [`DisjunctiveOracle`](@ref)
"""
mutable struct DisjunctiveRootNodePreprocessing <: AbstractRootNodePreprocessing
    typical_oracle::AbstractTypicalOracle
    disjunctive_oracle::AbstractDisjunctiveOracle
    seq_type::Type{<:AbstractBendersSeq}
    params::AbstractBendersSeqParam

    function DisjunctiveRootNodePreprocessing(
        typical_oracle::AbstractTypicalOracle,
        disjunctive_oracle::AbstractDisjunctiveOracle;
        seq_type::Type{<:AbstractBendersSeq} = BendersSeq,
        params::AbstractBendersSeqParam = BendersSeqParam()
    )
        new(typical_oracle, disjunctive_oracle, seq_type, params)
    end
end

"""
    root_node_processing!(master::AbstractMaster, preprocessing::DisjunctiveRootNodePreprocessing) -> Float64

Run a two-phase root-node preprocessing procedure using both the typical and
disjunctive oracles defined in a `DisjunctiveRootNodePreprocessing` object, and
return the total preprocessing time in seconds.

This routine performs the following steps:

1. **Model relaxation.**  
   All integrality constraints in `master.model` are temporarily relaxed.  
   The original integrality settings are restored automatically upon exit,
   even if an error occurs.

2. **Phase 1 – preprocessing with typical oracle.**  
   A preprocessing is built with `preprocessing.typical_oracle` and
   executed. 

4. **Phase 2 – preprocessing with disjunctive oracle.**  
   A preprocessing follows using `preprocessing.disjunctive_oracle`
   and solved using the remaining time budget.

The function returns the preprocessing times.
"""
function root_node_processing!(master::AbstractMaster, preprocessing::DisjunctiveRootNodePreprocessing)

    root_param = deepcopy(root_preprocessing.params)

    # Relax integrality, ensure undo() always runs even on error
    undo = relax_integrality(master.model)

    tic = time()
    try
        # Phase 1: Preprocessing with typical oracle
        env_root_typical = preprocessing.seq_type(master, preprocessing.typical_oracle; param = root_param)
        solve!(env_root_typical)
        typical_time = time() - tic

        root_param.time_limit -= typical_time

        # Phase 2: Preprocessing with disjunctive oracle 
        env_root_disjunctive = preprocessing.seq_type(data, master, preprocessing.disjunctive_oracle; param = root_param)
        solve!(env_root_disjunctive)
    finally
        # always restore integrality (even on exceptions)
        undo()
    end

    return time() - tic
end
