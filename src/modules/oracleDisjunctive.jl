export DisjunctiveOracle, DisjunctiveOracleParam

mutable struct DisjunctiveOracleParam <: AbstractOracleParam
    norm::AbstractNorm
    split_index_selection_rule::SplitIndexSelectionRule
    disjunctive_cut_append_rule::DisjunctiveCutsAppendRule
    strengthened::Bool
    add_benders_cuts_to_master::Bool
    fraction_of_benders_cuts_to_master::Float64
    reuse_dcglp::Bool
    lift::Bool 

    function DisjunctiveOracleParam(; 
                                    norm::AbstractNorm = LpNorm(Inf), 
                                    split_index_selection_rule::SplitIndexSelectionRule = RandomFractional(), disjunctive_cut_append_rule::DisjunctiveCutsAppendRule = AllDisjunctiveCuts(),
                                    strengthened::Bool = true, 
                                    add_benders_cuts_to_master::Bool = true, 
                                    fraction_of_benders_cuts_to_master::Float64 = 1.0, 
                                    reuse_dcglp::Bool = true,
                                    lift::Bool = false) 
        new(norm, split_index_selection_rule, disjunctive_cut_append_rule, strengthened, add_benders_cuts_to_master, fraction_of_benders_cuts_to_master, reuse_dcglp, lift)
    end
end 

function DisjunctiveOracle(data, 
                           typical_oracles::Vector{T}; 
                           param::DcglpParam = DcglpParam(),
                           solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9),
                           oracle_param::DisjunctiveOracleParam = DisjunctiveOracleParam()) where {T<:AbstractTypicalOracle}
    if typeof(oracle_param.norm) <: LpNorm
        return DisjunctiveOracleGammaNorm(data, typical_oracles, param, solver_param = solver_param, oracle_param = oracle_param)
    elseif typeof(oracle_param.norm) <: StandardNorm
        return DisjunctiveOracleStdNorm(data, typical_oracles, param, solver_param = solver_param, oracle_param = oracle_param)
    else
        throw(UndefError("DisjunctiveOracle: $(typeof(oracle_param.norm)) is not implemented."))
    end
end

"""
The following optional arguments are needed for `SpecializedBendersSeq`
- `throw_typical_cuts_for_errors` determines whether to return a typical Benders cut, when DCGLP encounters some issue. It must be `false` for SpecializedBendersSeq.
- `include_disjuctive_cuts_to_hyperplanes` determines whether to add a disjunctive cut found to `hyperplanes` to be returned; if it is `false`, the disjuctive cut should be added at a desired place via `oracle.disjunctiveCuts` or `oracle.disjunctiveCutsByIndex`. It must be `false` for SpecializedBendersSeq.
"""
function generate_cuts(oracle::AbstractDisjunctiveOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; rtol = 1e-6, atol = 1e-9, time_limit = 3600.0, throw_typical_cuts_for_errors = true, include_disjuctive_cuts_to_hyperplanes = true)

    tic = time()
    
    push!(oracle.splits, select_disjunctive_inequality(x_value, oracle.oracle_param.split_index_selection_rule))
    
    if get_sec_remaining(tic, time_limit) <= 0.0
        throw(TimeLimitException("Time limit reached during cut generation"))
    end

    replace_disjunctive_inequality!(oracle)
    
    # delete benders cuts previously added when not reusing dcglp
    if !oracle.oracle_param.reuse_dcglp
        if haskey(oracle.dcglp, :con_benders)
            delete.(oracle.dcglp, oracle.dcglp[:con_benders]) 
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

    # Retrieve zero and one indices if lifting is enabled
    zero_indices, one_indices = oracle.oracle_param.lift ? retrieve_zero_one(x_value, atol=atol) : (Int[], Int[]) 

    add_lifting_constraints!(oracle.dcglp, zero_indices, one_indices) 

    return solve_dcglp!(oracle, x_value, t_value, zero_indices, one_indices; rtol = rtol, atol = atol, time_limit = time_limit, throw_typical_cuts_for_errors = throw_typical_cuts_for_errors, include_disjuctive_cuts_to_hyperplanes = include_disjuctive_cuts_to_hyperplanes)
end

"""
Updates parameters of the DisjunctiveOracle. Changing the normalization updates the dcglp model, which is initially set during declaration.
"""
function set_parameter!(oracle::AbstractDisjunctiveOracle, param::DisjunctiveOracleParam)
    oracle.oracle_param = param
    if typeof(oracle) <: DisjunctiveOracleGammaNorm
        if haskey(oracle.dcglp, :concone)
            delete.(oracle.dcglp, oracle.dcglp[:concone]) 
            unregister(oracle.dcglp, :concone)
        end
        add_normalization_constraint(oracle.dcglp, oracle.oracle_param.norm)
    elseif typeof(oracle) <: DisjunctiveOracleStdNorm
        if typeof(param.norm) != StandardNorm
            @warn "set_parameter!: `StandardNorm` cannot be changed to `LpNorm`. If you intend to change it, please set it at the declaration of the oracle. Proceeding using `StandardNorm`..."
        end
    end
end

function set_parameter!(oracle::AbstractDisjunctiveOracle, param::String, value::Any)
    sym_param = Symbol(param)
    if sym_param ∈ fieldnames(typeof(oracle.oracle_param))
        setfield!(oracle.oracle_param, sym_param, value)
    else
        throw(ArgumentError("Parameter `$(param)` not found in `$(typeof(oracle.oracle_param))` for oracle of type `$(typeof(oracle))`"))
    end

    if sym_param == :norm
        if typeof(oracle) <: DisjunctiveOracleGammaNorm
            if haskey(oracle.dcglp, :concone)
                delete.(oracle.dcglp, oracle.dcglp[:concone]) 
                unregister(oracle.dcglp, :concone)
            end
            add_normalization_constraint(oracle.dcglp, oracle.oracle_param.norm)
        elseif typeof(oracle) <: DisjunctiveOracleStdNorm
            if typeof(param.norm) != StandardNorm
                @warn "set_parameter!: `StandardNorm` cannot be changed to `LpNorm`. If you intend to change it, please set it at the declaration of the oracle. Proceeding using `StandardNorm`..."
            end
        end
    end
end

"""
prototypes for user-customizable functions for DisjunctiveOracle
"""
function select_disjunctive_inequality(x_value::Vector{Float64}, split_selection_rule::SplitIndexSelectionRule; zero_tol = 1e-2)    
    throw(UndefError("update select_disjunctive_inequality for $(typeof(split_selection_rule))"))
    # should return a split: phi, phi_0
end

function add_disjunctive_cuts!(oracle::AbstractDisjunctiveOracle, rule::DisjunctiveCutsAppendRule)
    throw(UndefError("update add_disjunctive_cuts! for $(typeof(rule))"))
    # should add to dcglp
end

include("oracleDisjunctiveInterface.jl")

"""
utility functions for DisjunctiveOracle
"""
function get_split_index(oracle::AbstractDisjunctiveOracle)
    if !(typeof(oracle.oracle_param.split_index_selection_rule) <: SimpleSplit)
        throw(AlgorithmException("get_split_index is only valid for SimpleSplit"))
    end
    return findfirst(x -> x > 0.5, oracle.splits[end][1])
end

function replace_disjunctive_inequality!(oracle::AbstractDisjunctiveOracle)
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
    if typeof(oracle.oracle_param.norm) <: LpNorm
        @constraint(dcglp, con_split_kappa, 0 >= dcglp[:omega_0][1]*(phi_0+1) - phi' * dcglp[:omega_x][1,:])
        @constraint(dcglp, con_split_nu, 0 >= -dcglp[:omega_0][2]*phi_0 + phi' * dcglp[:omega_x][2,:])
    elseif typeof(oracle.oracle_param.norm) <: StandardNorm
        @constraint(dcglp, con_split_kappa, dcglp[:tau] >= dcglp[:omega_0][1]*(phi_0+1) - phi' * dcglp[:omega_x][1,:])
        @constraint(dcglp, con_split_nu, dcglp[:tau] >= -dcglp[:omega_0][2]*phi_0 + phi' * dcglp[:omega_x][2,:])
    else
        throw(UndefError("update replace_disjunctive_inequality! for $(typeof(oracle.oracle_param.norm))"))
    end
end

function solve_dcglp!(oracle::AbstractDisjunctiveOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; time_limit = time_limit)
    throw(UndefError("update solve_dcglp! for $(typeof(oracle))"))
end

function retrieve_zero_one(x_value::Vector{Float64}; atol = 1e-9)
    zeros_indices = findall(x -> isapprox(x, 0.0; atol=atol), x_value)
    ones_indices = findall(x -> isapprox(x, 1.0; atol=atol), x_value)
    return zeros_indices, ones_indices
end
# better to have oracle::AbstractDisjunctiveOracle, x_value; atol as inputs
# outputs zeros_indices, once_indices
function add_lifting_constraints!(dcglp::Model, zero_indices::Vector{Int}, one_indices::Vector{Int})
    # remove previously added lifting constraints
    haskey(dcglp, :con_zeta) && (delete.(dcglp, vcat(dcglp[:con_zeta]...)); unregister(dcglp, :con_zeta))
    haskey(dcglp, :con_xi) && (delete.(dcglp, vcat(dcglp[:con_xi]...)); unregister(dcglp, :con_xi))

    # add lifting constraints
    !isempty(zero_indices) && @constraint(dcglp, con_zeta[i in 1:2, j=1:length(zero_indices)], 0 >= dcglp[:omega_x][i, zero_indices[j]])
    !isempty(one_indices) && @constraint(dcglp, con_xi[i in 1:2, j=1:length(one_indices)], 0 >= dcglp[:omega_0][i] - dcglp[:omega_x][i, one_indices[j]])
end