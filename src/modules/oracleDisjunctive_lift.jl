export DisjunctiveOracle, DisjunctiveOracleParam

mutable struct DisjunctiveOracleParam <: AbstractOracleParam
    norm::AbstractNorm
    split_index_selection_rule::SplitIndexSelectionRule
    disjunctive_cut_append_rule::DisjunctiveCutsAppendRule
    strengthened::Bool
    add_benders_cuts_to_master::Bool
    fraction_of_benders_cuts_to_master::Float64
    reuse_dcglp::Bool
    lift::Bool # lifting

    function DisjunctiveOracleParam(; 
                                    norm::AbstractNorm = LpNorm(Inf), 
                                    split_index_selection_rule::SplitIndexSelectionRule = RandomFractional(), disjunctive_cut_append_rule::DisjunctiveCutsAppendRule = AllDisjunctiveCuts(),
                                    strengthened::Bool=true, 
                                    add_benders_cuts_to_master::Bool=true, 
                                    fraction_of_benders_cuts_to_master::Float64 = 1.0, 
                                    reuse_dcglp::Bool=true,
                                    lift::Bool=true) # lifting
        new(norm, split_index_selection_rule, disjunctive_cut_append_rule, strengthened, add_benders_cuts_to_master, fraction_of_benders_cuts_to_master, reuse_dcglp, lift) # lifting
    end
end 

mutable struct DisjunctiveOracle <: AbstractDisjunctiveOracle
    
    oracle_param::DisjunctiveOracleParam

    dcglp::Model
    typical_oracles::Vector{AbstractTypicalOracle}

    # Dcglp loop parameters
    param::DcglpParam

    # log for splits and disjunctive cuts
    disjunctiveCutsByIndex::Vector{Vector{Hyperplane}}
    disjunctiveCuts::Vector{Hyperplane}
    splits::Vector{Tuple{SparseVector{Float64, Int}, Float64}}

    function DisjunctiveOracle(data, 
                               typical_oracles::Vector{T}; 
                               param::DcglpParam = DcglpParam(),
                               solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9),
                               oracle_param::DisjunctiveOracleParam = DisjunctiveOracleParam()) where {T<:AbstractTypicalOracle}
        @debug "Building disjunctive oracle"
        dcglp = Model()
        # Define variables
        @variable(dcglp, tau)
        @variable(dcglp, omega_0[1:2]) # 1 for kappa; 2 for nu
        @variable(dcglp, omega_x[1:2,1:data.dim_x])
        @variable(dcglp, omega_t[1:2,1:data.dim_t])
        @variable(dcglp, sx[1:data.dim_x])
        @variable(dcglp, st[1:data.dim_t])
        
        # Set objective
        @objective(dcglp, Min, tau)

        # Add constraints
        @constraint(dcglp, [i=1:2], omega_t[i,:] .>= -1e6 * omega_0[i])
        @constraint(dcglp, coneta[i in 1:2, j in 1:data.dim_x], 0 >= -omega_0[i] + omega_x[i,j]) 
        @constraint(dcglp, condelta[i in 1:2, j in 1:data.dim_x], 0 >= -omega_x[i,j])
        @constraint(dcglp, conineq[i in 1:2], omega_0[i] >= 0)

        # Add gamma constraints
        @constraint(dcglp, con0, omega_0[1] + omega_0[2] == 1)
        @constraint(dcglp, conx, omega_x[1,:] + omega_x[2,:] - sx .== 0)
        @constraint(dcglp, cont[j=1:data.dim_t], omega_t[1,j] + omega_t[2,j] - st[j] == 0) # must be in this form to recognize it as a vector

        assign_attributes!(dcglp, solver_param)
        
        add_normalization_constraint(dcglp, oracle_param.norm)
        
        disjunctiveCutsByIndex = [Vector{Hyperplane}() for i=1:data.dim_x]
        splits = Vector{Tuple{SparseVector{Float64, Int}, Float64}}()

        new(oracle_param, dcglp, typical_oracles, param, disjunctiveCutsByIndex, Vector{Hyperplane}(), splits)
    end
end 

function generate_cuts(oracle::DisjunctiveOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol = 1e-6, time_limit = 3600)
    _, _, t_value = generate_cuts(oracle.typical_oracles[1], x_value, t_value)
    tic = time()

    # Retrieve zero and one indices if lifting is enabled
    zero_indices, one_indices = oracle.oracle_param.lift ? retrieve_zero_one(x_value) : (Int[], Int[]) # lifting
    
    push!(oracle.splits, select_disjunctive_inequality(x_value, oracle.oracle_param.split_index_selection_rule))
    
    if get_sec_remaining(tic, time_limit) <= 0.0
        throw(TimeLimitException("Time limit reached during cut generation"))
    end

    replace_disjunctive_inequality!(oracle)
    
    # delete benders cuts previously added when not reusing dcglp
    if !oracle.oracle_param.reuse_dcglp
        if haskey(oracle.dcglp, :con_benders)
            delete(oracle.dcglp, oracle.dcglp[:con_benders]) 
            unregister(oracle.dcglp, :con_benders)
        end
    end

    # add previously found disjunctive cuts based on a user-given append rule
    add_disjunctive_cuts!(oracle, oracle.oracle_param.disjunctive_cut_append_rule)

    if get_sec_remaining(tic, time_limit) <= 0.0
        throw(TimeLimitException("Time limit reached during cut generation"))
    end

    set_normalized_rhs.(oracle.dcglp[:conx], x_value)
    set_normalized_rhs.(oracle.dcglp[:cont], t_value)

    @debug "zero_indices", zero_indices
    @debug "one_indices", one_indices

    add_lifting_constraints!(oracle.dcglp, zero_indices, one_indices) # lifting

    return solve_dcglp!(oracle, x_value, t_value, zero_indices, one_indices; time_limit = time_limit)
end
"""
Updates parameters of the DisjunctiveOracle. Changing the normalization updates the dcglp model, which is initially set during declaration.
"""
function set_parameter!(oracle::DisjunctiveOracle, param::DisjunctiveOracleParam)
    oracle.oracle_param = param
    if haskey(oracle.dcglp, :concone)
        delete(oracle.dcglp, oracle.dcglp[:concone]) 
        unregister(oracle.dcglp, :concone)
    end
    add_normalization_constraint(oracle.dcglp, oracle.oracle_param.norm)
end
  
function set_parameter!(oracle::DisjunctiveOracle, param::String, value::Any)
    sym_param = Symbol(param)
    if sym_param ∈ fieldnames(typeof(oracle.oracle_param))
        setfield!(oracle.oracle_param, sym_param, value)
    else
        throw(ArgumentError("Parameter `$(param)` not found in `$(typeof(oracle.oracle_param))` for oracle of type `$(typeof(oracle))`"))
    end

    if sym_param == :norm
        if haskey(oracle.dcglp, :concone)
            delete(oracle.dcglp, oracle.dcglp[:concone]) 
            unregister(oracle.dcglp, :concone)
        end
        add_normalization_constraint(oracle.dcglp, oracle.oracle_param.norm)
    end
end

"""
prototypes for user-customizable functions for DisjunctiveOracle
"""
function add_normalization_constraint(dcglp::Model, norm::AbstractNorm)
    throw(UndefError("update add_normalization_constraint for $(typeof(norm))"))
    # should add a normalization constraint to dcglp. 
end

function select_disjunctive_inequality(x_value::Vector{Float64}, split_selection_rule::SplitIndexSelectionRule; zero_tol = 1e-2)    
    throw(UndefError("update select_disjunctive_inequality for $(typeof(split_selection_rule))"))
    # should return a split: phi, phi_0
end

function add_disjunctive_cuts!(oracle::DisjunctiveOracle, rule::DisjunctiveCutsAppendRule)
    throw(UndefError("update add_disjunctive_cuts! for $(typeof(rule))"))
    # should add to dcglp
end

include("oracleDisjunctiveInterface.jl")

"""
utility functions for DisjunctiveOracle
"""
function get_split_index(oracle::DisjunctiveOracle)
    if !(typeof(oracle.oracle_param.split_index_selection_rule) <: SimpleSplit)
        throw(AlgorithmException("get_split_index is only valid for SimpleSplit"))
    end
    return findfirst(x -> x > 0.5, oracle.splits[end][1])
end

function replace_disjunctive_inequality!(oracle::DisjunctiveOracle)
    dcglp = oracle.dcglp
    phi = oracle.splits[end][1]
    phi_0 = oracle.splits[end][2]
    
    if haskey(dcglp, :con_split_kappa)
        delete(dcglp, dcglp[:con_split_kappa]) 
        unregister(dcglp, :con_split_kappa)
    end
    if haskey(dcglp, :con_split_nu)
        delete(dcglp, dcglp[:con_split_nu]) 
        unregister(dcglp, :con_split_nu)
    end
        
    # Add new constraints
    @constraint(dcglp, con_split_kappa, 0 >= dcglp[:omega_0][1]*(phi_0+1) - phi' * dcglp[:omega_x][1,:])
    @constraint(dcglp, con_split_nu, 0 >= -dcglp[:omega_0][2]*phi_0 + phi' * dcglp[:omega_x][2,:])
end

function solve_dcglp!(oracle::AbstractDisjunctiveOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; time_limit = time_limit)
    throw(UndefError("update solve_dcglp! for $(typeof(oracle))"))
end

function retrieve_zero_one(x_value::Vector{Float64}) # lifting
    # atol can be assigend from oracle_param
    zeros_indices = findall(x -> isapprox(x, 0.0; atol=1e-6), x_value)
    ones_indices = findall(x -> isapprox(x, 1.0; atol=1e-6), x_value)
    return zeros_indices, ones_indices
end

function add_lifting_constraints!(dcglp::Model, zero_indices::Vector{Int}, one_indices::Vector{Int}) # lifting
    # remove previously added lifting constraints
    haskey(dcglp, :con_zeta) && (delete.(dcglp, vcat(dcglp[:con_zeta]...)); unregister(dcglp, :con_zeta))
    haskey(dcglp, :con_xi) && (delete.(dcglp, vcat(dcglp[:con_xi]...)); unregister(dcglp, :con_xi))

    # add lifting constraints

    # dcglp[:con_zeta] = Vector{ConstraintRef}[ConstraintRef[] for _ in 1:2]
    # dcglp[:con_xi] = Vector{ConstraintRef}[ConstraintRef[] for _ in 1:2]
    # println("dual of empty", dual.(dcglp[:con_zeta][1,:]))

    !isempty(zero_indices) && @constraint(dcglp, con_zeta[i in 1:2, j=1:length(zero_indices)], 0>= dcglp[:omega_x][i, zero_indices[j]])
    !isempty(one_indices) && @constraint(dcglp, con_xi[i in 1:2, j=1:length(one_indices)], 0>= dcglp[:omega_0][i] - dcglp[:omega_x][i, one_indices[j]])
end