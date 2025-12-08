export DisjunctiveOracle, DisjunctiveOracleParam

"""
    DisjunctiveOracleParam <: AbstractOracleParam

Parameters for configuring a DisjunctiveOracle in Benders decomposition.

This structure contains all the parameters needed to control the behavior of the disjunctive cut generation process via the DCGLP (Disjunctive Cut Generating Linear Program).

# Fields
- `norm::AbstractNorm`: Norm type used for normalization in DCGLP (default: `LpNorm(Inf)`)
- `split_index_selection_rule::SplitIndexSelectionRule`: Rule for selecting which variable to split on (default: `RandomFractional()`)
- `disjunctive_cut_append_rule::DisjunctiveCutsAppendRule`: Rule for adding previously found disjunctive cuts to DCGLP (default: `AllDisjunctiveCuts()`)
- `strengthened::Bool`: Whether to apply strengthening to the disjunctive cuts (default: `true`)
- `add_benders_cuts_to_master::Int`: Controls adding Benders cuts to master problem (default: `1`)
  - `0`: Do not add Benders cuts to master
  - `1`: Add all Benders cuts regardless of violation at (x,t)
  - `2`: Add only violated Benders cuts
- `fraction_of_benders_cuts_to_master::Float64`: Fraction of Benders cuts to add to master (default: `1.0`)
- `reuse_dcglp::Bool`: Whether to reuse DCGLP model across iterations (default: `true`)
- `lift::Bool`: Whether to use lifting for variables fixed to 0 or 1 (default: `false`)
- `adjust_t_to_fx::Bool`: Whether to adjust t to f(x) before solving DCGLP (default: `false`)
- `zero_tol::Float64`: Tolerance for considering values as zero (default: `1e-9`)

# Examples
```julia
# Create with default parameters
param = DisjunctiveOracleParam()

# Create with custom parameters
param = DisjunctiveOracleParam(
    norm = LpNorm(1.0),
    split_index_selection_rule = MostFractional(),
    strengthened = true,
    add_benders_cuts_to_master = 2  # Only violated cuts
)
```

See also: [`DisjunctiveOracle`](@ref), [`LpNorm`](@ref)
"""
mutable struct DisjunctiveOracleParam <: AbstractOracleParam
    
    dcglp_param::DcglpParam # Dcglp loop parameters
    norm::AbstractNorm
    split_index_selection_rule::SplitIndexSelectionRule
    disjunctive_cut_append_rule::DisjunctiveCutsAppendRule
    strengthened::Bool
    add_benders_cuts_to_master::Int # 0: do not add; 1: add regardless of violation at (x,t); 2: only violated ones
    fraction_of_benders_cuts_to_master::Float64
    reuse_dcglp::Bool
    lift::Bool 
    adjust_t_to_fx::Bool
    zero_tol::Float64
    

    function DisjunctiveOracleParam(dcglp_param::DcglpParam; 
                                    norm::AbstractNorm = LpNorm(Inf), 
                                    split_index_selection_rule::SplitIndexSelectionRule = RandomFractional(), disjunctive_cut_append_rule::DisjunctiveCutsAppendRule = AllDisjunctiveCuts(),
                                    strengthened::Bool=true, 
                                    add_benders_cuts_to_master::Union{Bool,Int} = 1, 
                                    fraction_of_benders_cuts_to_master::Float64 = 1.0, 
                                    reuse_dcglp::Bool=true,
                                    lift::Bool=false,
                                    adjust_t_to_fx::Bool=false,
                                    zero_tol=1e-9) 
        add_bcuts_to_master = add_benders_cuts_to_master === true ? 1 : add_benders_cuts_to_master === false ? 0 : add_benders_cuts_to_master in (0, 1, 2) ? add_benders_cuts_to_master : throw(ArgumentError("`add_benders_cuts_to_master` must be true, false, or an integer in {0, 1, 2}"))
        
        new(dcglp_param, norm, split_index_selection_rule, disjunctive_cut_append_rule, strengthened, add_bcuts_to_master, fraction_of_benders_cuts_to_master, reuse_dcglp, lift, adjust_t_to_fx, zero_tol)
    end
end

"""
    DisjunctiveOracle <: AbstractDisjunctiveOracle

An oracle for generating disjunctive cuts in Benders decomposition using the DCGLP method.

The DisjunctiveOracle combines two typical oracles (for the kappa and nu subproblems) and
solves a Disjunctive Cut Generating Linear Program (DCGLP) to produce stronger cuts than
standard Benders cuts. This can significantly improve convergence for mixed-integer programs.

# Fields
- `oracle_param::DisjunctiveOracleParam`: Parameters controlling the oracle's behavior
- `dcglp::Model`: The JuMP model for the Disjunctive Cut Generating Linear Program
- `typical_oracles::Vector{AbstractTypicalOracle}`: Two typical oracles (index 1: kappa, index 2: nu)
- `param::DcglpParam`: Parameters for the DCGLP cutting-plane loop
- `disjunctiveCutsByIndex::Vector{Vector{Hyperplane}}`: Disjunctive cuts organized by split variable index
- `disjunctiveCuts::Vector{Hyperplane}`: All disjunctive cuts generated so far
- `splits::Vector{Tuple{SparseVector{Float64, Int}, Float64}}`: History of split inequalities (phi, phi_0)

# Constructor
```julia
DisjunctiveOracle(data, typical_oracles::Vector{T};
                  param::DcglpParam = DcglpParam(),
                  solver_param::Dict{String,Any} = Dict(...),
                  oracle_param::DisjunctiveOracleParam = DisjunctiveOracleParam()) where {T<:AbstractTypicalOracle}
```

# Arguments
- `data`: Problem data containing dimensions
- `typical_oracles`: Vector of two typical oracles for kappa and nu subproblems
- `param`: DCGLP loop parameters (optional)
- `solver_param`: Solver configuration for DCGLP model (optional)
- `oracle_param`: Oracle behavior parameters (optional)

# Examples
```julia
# Create typical oracles
oracle_kappa = ClassicalOracle(data, scen_idx=1)
oracle_nu = ClassicalOracle(data, scen_idx=2)

# Create disjunctive oracle with default parameters
disj_oracle = DisjunctiveOracle(data, [oracle_kappa, oracle_nu])

# Create with custom parameters
disj_param = DisjunctiveOracleParam(strengthened = true, lift = true)
disj_oracle = DisjunctiveOracle(data, [oracle_kappa, oracle_nu]; oracle_param = disj_param)
```

# Notes
- The DCGLP model is built automatically during construction
- The oracle maintains a history of all disjunctive cuts and splits
- Disjunctive cuts can be organized by split variable index for specialized algorithms

See also: [`DisjunctiveOracleParam`](@ref), [`generate_cuts`](@ref), [`solve_dcglp!`](@ref)
"""
mutable struct DisjunctiveOracle <: AbstractDisjunctiveOracle
    
    oracle_param::DisjunctiveOracleParam

    dcglp::Model
    typical_oracles::Vector{AbstractTypicalOracle}

    # log for splits and disjunctive cuts
    disjunctiveCutsByIndex::Vector{Vector{Hyperplane}}
    disjunctiveCuts::Vector{Hyperplane}
    splits::Vector{Tuple{SparseVector{Float64, Int}, Float64}}

    function DisjunctiveOracle(data, 
                               typical_oracles::Vector{T}; 
                               param::DisjunctiveOracleParam = DisjunctiveOracleParam()) where {T<:AbstractTypicalOracle}
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

    # oracle_param should not be optional unless we have default software-free optimizer
    function DisjunctiveOracle(master::Master, 
                            typical_oracles::Vector{T},
                            oracle_param::DisjunctiveOracleParam) where {T<:AbstractTypicalOracle}
        @debug "Building disjunctive oracle"

        # Initialize data object
        dim_x = length(master.x)
        obj = objective_function(master.model)
        c_x = [coefficient(obj, master.x[i]) for i in 1:dim_x]
        dim_t = length(master.t)
        c_t = [coefficient(obj, master.t[i]) for i in 1:dim_t]
        data = Data(dim_x, dim_t, EmptyData(), c_x, c_t)

        for xi in master.x
            if !is_binary(xi)
                @error "Split oracles currently require all master variables to be binary."
            end
        end

        # Initialize dcglp problem
        dcglp = Model(oracle_param.dcglp_param.optimizer)
        
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

        for i=1:2
            transfer_scaled_linear_rows_and_bounds_with_types!(master.model, master.x, dcglp, omega_x[i,:], omega_0[i])
        end

        add_normalization_constraint(dcglp, oracle_param.norm)

        disjunctiveCutsByIndex = [Vector{Hyperplane}() for i=1:data.dim_x]
        splits = Vector{Tuple{SparseVector{Float64, Int}, Float64}}()

        new(oracle_param, dcglp, typical_oracles, disjunctiveCutsByIndex, Vector{Hyperplane}(), splits)
    end
end

"""
    generate_cuts(oracle::DisjunctiveOracle, x_value, t_value; kwargs...) -> (Bool, Vector{Hyperplane}, Vector{Float64})

Generate cuts for a disjunctive oracle by solving the DCGLP problem.

This function selects a disjunctive inequality, updates the DCGLP model, and calls `solve_dcglp!` to generate either a disjunctive cut or fall back to typical Benders cuts.

# Arguments
- `oracle::DisjunctiveOracle`: The disjunctive oracle containing DCGLP model and typical oracles
- `x_value::Vector{Float64}`: Current first-stage solution
- `t_value::Vector{Float64}`: Current second-stage approximation

# Keyword Arguments
- `time_limit::Float64`: Maximum time allowed for cut generation (default: 3600.0)
- `throw_typical_cuts_for_errors::Bool`: If true, return typical Benders cuts when DCGLP encounters errors; if false, throw an exception instead (default: true)
- `include_disjunctive_cuts_to_hyperplanes::Bool`: If true, add the generated disjunctive cut to the returned hyperplanes; if false, the cut is only stored in `oracle.disjunctiveCuts` (default: true)

# Returns
A tuple `(is_in_L, hyperplanes, f_x)`:
- `is_in_L::Bool`: Whether the point is feasible (false if cuts were generated)
- `hyperplanes::Vector{Hyperplane}`: Generated cuts to add to the master problem
- `f_x::Vector{Float64}`: Subproblem objective values

# Notes
- The parameter `throw_typical_cuts_for_errors` must be set to `false` when using `SpecializedBendersSeq`
- When `include_disjunctive_cuts_to_hyperplanes` is `false`, the disjunctive cut can be accessed via `oracle.disjunctiveCuts` or `oracle.disjunctiveCutsByIndex` for specialized algorithms

See also: [`solve_dcglp!`](@ref), [`DisjunctiveOracle`](@ref)
"""
function generate_cuts(oracle::DisjunctiveOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; time_limit = 3600.0, throw_typical_cuts_for_errors = true, include_disjunctive_cuts_to_hyperplanes = true)

    tic = time()
    
    push!(oracle.splits, select_disjunctive_inequality(x_value, oracle.oracle_param.split_index_selection_rule; zero_tol = oracle.oracle_param.zero_tol))
    
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
    zero_indices, one_indices = oracle.oracle_param.lift ? retrieve_zero_one(x_value, oracle.oracle_param.zero_tol) : (Int[], Int[])

    add_lifting_constraints!(oracle.dcglp, zero_indices, one_indices) 

    return solve_dcglp!(oracle, x_value, t_value, zero_indices, one_indices; time_limit = time_limit, throw_typical_cuts_for_errors = throw_typical_cuts_for_errors, include_disjunctive_cuts_to_hyperplanes = include_disjunctive_cuts_to_hyperplanes)
end
"""
Updates parameters of the DisjunctiveOracle. Changing the normalization updates the dcglp model, which is initially set during declaration.
"""
function set_parameter!(oracle::DisjunctiveOracle, param::DisjunctiveOracleParam)
    oracle.oracle_param = param
    if haskey(oracle.dcglp, :concone)
        delete.(oracle.dcglp, oracle.dcglp[:concone]) 
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
            delete.(oracle.dcglp, oracle.dcglp[:concone]) 
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
include("Dcglp.jl")

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

function retrieve_zero_one(x_value::Vector{Float64}, zero_tol)
    zeros_indices = findall(x -> isapprox(x, 0.0; atol=zero_tol), x_value)
    ones_indices = findall(x -> isapprox(x, 1.0; atol=zero_tol), x_value)
    return zeros_indices, ones_indices
end

function add_lifting_constraints!(dcglp::Model, zero_indices::Vector{Int}, one_indices::Vector{Int})
    # remove previously added lifting constraints
    haskey(dcglp, :con_zeta) && (delete.(dcglp, vcat(dcglp[:con_zeta]...)); unregister(dcglp, :con_zeta))
    haskey(dcglp, :con_xi) && (delete.(dcglp, vcat(dcglp[:con_xi]...)); unregister(dcglp, :con_xi))

    # add lifting constraints
    !isempty(zero_indices) && @constraint(dcglp, con_zeta[i in 1:2, j=1:length(zero_indices)], 0>= dcglp[:omega_x][i, zero_indices[j]])
    !isempty(one_indices) && @constraint(dcglp, con_xi[i in 1:2, j=1:length(one_indices)], 0>= dcglp[:omega_0][i] - dcglp[:omega_x][i, one_indices[j]])
end