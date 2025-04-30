export DisjunctiveOracleGammaNorm

mutable struct DisjunctiveOracleGammaNorm <: AbstractDisjunctiveOracle
    
    oracle_param::DisjunctiveOracleParam

    dcglp::Model
    typical_oracles::Vector{AbstractTypicalOracle}

    # Dcglp loop parameters
    param::DcglpParam

    # log for splits and disjunctive cuts
    disjunctiveCutsByIndex::Vector{Vector{Hyperplane}}
    disjunctiveCuts::Vector{Hyperplane}
    splits::Vector{Tuple{SparseVector{Float64, Int}, Float64}}

    function DisjunctiveOracleGammaNorm(data, 
                                      typical_oracles::Vector{T}; 
                                      param::DcglpParam = DcglpParam(),
                                      solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9),
                                      oracle_param::DisjunctiveOracleParam = DisjunctiveOracleParam()) where {T<:AbstractTypicalOracle}
        @debug "Building disjunctive oracle"
        dcglp = Model()
        # Define variables
        @variable(dcglp, tau)
        @variable(dcglp, omega_0[1:2]) # 1 for kappa; 2 for nu
        @variable(dcglp, omega_x[1:2,1:dim_x])
        @variable(dcglp, omega_t[1:2,1:dim_t])
        @variable(dcglp, sx[1:dim_x])
        @variable(dcglp, st[1:dim_t])
        
        # Set objective
        @objective(dcglp, Min, tau)
    
        # Add constraints
        @constraint(dcglp, [i=1:2], omega_t[i,:] .>= -1e6 * omega_0[i])
        @constraint(dcglp, coneta[i in 1:2, j in 1:dim_x], 0 >= -omega_0[i] + omega_x[i,j]) 
        @constraint(dcglp, condelta[i in 1:2, j in 1:dim_x], 0 >= -omega_x[i,j])
        @constraint(dcglp, conineq[i in 1:2], omega_0[i] >= 0)
    
        # Add gamma constraints
        @constraint(dcglp, con0, omega_0[1] + omega_0[2] == 1)
        @constraint(dcglp, conx, omega_x[1,:] + omega_x[2,:] - sx .== 0)
        @constraint(dcglp, cont[j=1:dim_t], omega_t[1,j] + omega_t[2,j] - st[j] == 0) # must be in this form to recognize it as a vector
        
        add_normalization_constraint(dcglp, norm)

        assign_attributes!(dcglp, solver_param)
        
        disjunctiveCutsByIndex = [Vector{Hyperplane}() for i=1:data.dim_x]
        splits = Vector{Tuple{SparseVector{Float64, Int}, Float64}}()

        new(oracle_param, dcglp, typical_oracles, param, disjunctiveCutsByIndex, Vector{Hyperplane}(), splits)
    end
end

"""
utility functions for DisjunctiveOracleGammaNorm
"""
function add_normalization_constraint(dcglp::Model, norm::LpNorm)
    # CPLEX only accepts p=1,2,Inf
    # if conic solver, we can use the following line
    # @constraint(dcglp, concone, var_vec in MOI.NormCone(norm.p, data.dim_x + data.dim_t + 1))
    var_vec = [dcglp[:tau]; dcglp[:sx]; dcglp[:st]]
    
    if norm.p == 1.0
        @constraint(dcglp, concone, var_vec in MOI.NormOneCone(length(var_vec)))
    elseif norm.p == 2.0
        @constraint(dcglp, concone, var_vec in MOI.SecondOrderCone(length(var_vec)))
    elseif norm.p == Inf
        @constraint(dcglp, concone, var_vec in MOI.NormInfinityCone(length(var_vec)))
    else
        throw(UndefError("Unsupported LpNorm: p=$(norm.p)"))
    end
end