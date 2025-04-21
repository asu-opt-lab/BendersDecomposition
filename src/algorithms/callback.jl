export RootNodePreprocessing, LazyCallback, UserCallback
export lazy_callback, user_callback
export EmptyCallbackParam, UserCallbackParam
export NoRootNodePreprocessing, NoUserCallback

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

"""
    AbstractCallbackParam

Abstract type for parameters used in callbacks during the branch-and-bound process.
"""
abstract type AbstractCallbackParam end

"""
    EmptyCallbackParam <: AbstractCallbackParam

Represents empty (default) parameters for callbacks.
"""
struct EmptyCallbackParam <: AbstractCallbackParam
end

"""
    AbstractLazyCallback

Abstract type for lazy constraint callbacks in Benders decomposition.
"""
abstract type AbstractLazyCallback end

"""
    LazyCallback <: AbstractLazyCallback

Configuration for lazy constraint callbacks in the branch-and-bound process.
Used to dynamically add Benders cuts when integer solutions are found.

# Fields
- `params::AbstractCallbackParam`: Parameters for the callback
- `oracle::AbstractOracle`: Oracle used to generate Benders cuts
"""
struct LazyCallback <: AbstractLazyCallback
    params::AbstractCallbackParam
    oracle::AbstractOracle
    

    function LazyCallback(; params=EmptyCallbackParam(), oracle)
        new(params, oracle)
    end
    
    function LazyCallback(data; params=EmptyCallbackParam())
        new(params, ClassicalOracle(data))
    end
end

"""
    UserCallbackParam <: AbstractCallbackParam

Parameters for user callbacks in the branch-and-bound process.

# Fields
- `frequency::Int = 50`: How often to process nodes (every N fractional nodes)
- `node_count::Int = -1`: Only process nodes after this node count (-1 means process all)
- `depth::Int = -1`: Only process nodes with depth >= this value (-1 means process all depths)
"""
Base.@kwdef struct UserCallbackParam <: AbstractCallbackParam
    frequency::Int = 50
    node_count::Int = -1
    depth::Int = -1
end

"""
    AbstractUserCallback

Abstract type for user cut callbacks in Benders decomposition.
"""
abstract type AbstractUserCallback end

"""
    NoUserCallback <: AbstractUserCallback

Indicates that no user callbacks should be used.
"""
struct NoUserCallback <: AbstractUserCallback end

"""
    UserCallback <: AbstractUserCallback

Configuration for user cut callbacks in the branch-and-bound process.
Used to dynamically add Benders cuts at fractional nodes.

# Fields
- `params::UserCallbackParam`: Parameters controlling when cuts are generated
- `oracle::AbstractOracle`: Oracle used to generate Benders cuts
"""
struct UserCallback <: AbstractUserCallback
    params::UserCallbackParam
    oracle::AbstractOracle
    
    function UserCallback(; params=UserCallbackParam(), oracle)
        new(params, oracle)
    end
    
    function UserCallback(data; params=UserCallbackParam())
        new(params, ClassicalOracle(data))
    end
end

"""
    lazy_callback(cb_data, master_model::Model, log::BendersBnBLog, callback::LazyCallback)

Callback function for adding lazy constraints in the branch-and-bound process.
Generates and adds Benders cuts when integer solutions are found.

# Arguments
- `cb_data`: Callback data from the solver
- `master_model::Model`: The JuMP master problem model
- `log::BendersBnBLog`: Log object to record statistics
- `callback::LazyCallback`: Configuration for the lazy callback
"""
function lazy_callback(cb_data, master_model::Model, log::BendersBnBLog, callback::LazyCallback)
    status = JuMP.callback_node_status(cb_data, master_model)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER

        state = BendersBnBState()
        if solver_name(master_model) == "CPLEX"
            n_count = Ref{CPXINT}()
            ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)
            state.node = n_count[]
        end

        state.values[:x] = JuMP.callback_value.(cb_data, master_model[:x])
        state.values[:t] = JuMP.callback_value.(cb_data, master_model[:t])

        state.oracle_time = @elapsed begin
            state.is_in_L, hyperplanes, state.f_x = generate_cuts(callback.oracle, state.values[:x], state.values[:t])
            cuts = !state.is_in_L ? hyperplanes_to_expression(master_model, hyperplanes, master_model[:x], master_model[:t]) : []
        end

        if !isempty(cuts)
            for cut in cuts
                cut_constraint = @build_constraint(0 >= cut)
                MOI.submit(master_model, MOI.LazyConstraint(cb_data), cut_constraint)
                state.num_cuts += 1
            end
        end
        record_node!(log, state, true)
    end
end

"""
    user_callback(cb_data, master_model::Model, log::BendersBnBLog, callback::UserCallback)

Callback function for adding user cuts in the branch-and-bound process.
Generates and adds Benders cuts at fractional nodes based on the specified frequency and criteria.

# Arguments
- `cb_data`: Callback data from the solver
- `master_model::Model`: The JuMP master problem model
- `log::BendersBnBLog`: Log object to record statistics
- `callback::UserCallback`: Configuration for the user callback with parameters controlling when cuts are generated
"""
function user_callback(cb_data, master_model::Model, log::BendersBnBLog, callback::UserCallback)
    status = JuMP.callback_node_status(cb_data, master_model)
    
    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        log.num_of_fraction_node += 1
        
        # Check if we should process this node based on frequency
        if log.num_of_fraction_node >= callback.params.frequency
            log.num_of_fraction_node = 0
            
            # Get node information if using CPLEX
            process_node = true
            if solver_name(master_model) == "CPLEX" && (callback.params.node_count != -1 || callback.params.depth != -1)
                n_count = Ref{CPXINT}()
                node_depth = Ref{CPXINT}()
                CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)
                CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODEDEPTH, node_depth)
                
                # Check if node meets criteria
                if (callback.params.node_count != -1 && n_count[] > callback.params.node_count) || 
                   (callback.params.depth != -1 && node_depth[] < callback.params.depth)
                    process_node = false
                end
            elseif (callback.params.node_count != -1 || callback.params.depth != -1) && solver_name(master_model) != "CPLEX"
                @warn "node_count and depth parameters are not supported for $(solver_name(master_model)) solver"
            end
            
            if process_node
                # Create state and get current variable values
                state = BendersBnBState()
                state.values[:x] = JuMP.callback_value.(cb_data, master_model[:x])
                state.values[:t] = JuMP.callback_value.(cb_data, master_model[:t])
                
                # Generate cuts
                state.oracle_time = @elapsed begin
                    state.is_in_L, hyperplanes, state.f_x = generate_cuts(callback.oracle, state.values[:x], state.values[:t])
                    cuts = !state.is_in_L ? hyperplanes_to_expression(master_model, hyperplanes, master_model[:x], master_model[:t]) : []
                end
                
                # Add cuts if any were generated
                if !isempty(cuts)
                    for cut in cuts
                        cut_constraint = @build_constraint(0 >= cut)
                        MOI.submit(master_model, MOI.UserCut(cb_data), cut_constraint)
                        state.num_cuts += 1
                    end
                end
                
                # Record node information
                record_node!(log, state, false)
            end
        end
    end
end