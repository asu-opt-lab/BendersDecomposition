export DisjunctiveOracleStdNorm

mutable struct DisjunctiveOracleStdNorm <: AbstractDisjunctiveOracle
    
    oracle_param::DisjunctiveOracleParam

    dcglp::Model
    typical_oracles::Vector{AbstractTypicalOracle}

    # Dcglp loop parameters
    param::DcglpParam

    # log for splits and disjunctive cuts
    disjunctiveCutsByIndex::Vector{Vector{Hyperplane}}
    disjunctiveCuts::Vector{Hyperplane}
    splits::Vector{Tuple{SparseVector{Float64, Int}, Float64}}

    function DisjunctiveOracleStdNorm(data, 
                                      typical_oracles::Vector{T}; 
                                      param::DcglpParam = DcglpParam(),
                                      solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9),
                                      oracle_param::DisjunctiveOracleParam = DisjunctiveOracleParam()) where {T<:AbstractTypicalOracle}
        @debug "Building disjunctive oracle"
        dcglp = Model()
        # Define variables
        @variable(dcglp, tau >= 0)
        @variable(dcglp, omega_0[1:2]) # 1 for kappa; 2 for nu
        @variable(dcglp, omega_x[1:2,1:dim_x])
        @variable(dcglp, omega_t[1:2,1:dim_t])
        
        # Set objective
        @objective(dcglp, Min, tau)
    
        # Add constraints
        @constraint(dcglp, [i=1:2], omega_t[i,:] .+ tau .>= -1e6 * omega_0[i])
        @constraint(dcglp, coneta[i in 1:2, j in 1:dim_x], tau >= -omega_0[i] + omega_x[i,j]) 
        @constraint(dcglp, condelta[i in 1:2, j in 1:dim_x], tau >= -omega_x[i,j])
        @constraint(dcglp, conineq[i in 1:2], omega_0[i] >= 0)
    
        # Add gamma constraints
        @constraint(dcglp, con0, omega_0[1] + omega_0[2] == 1)
        @constraint(dcglp, conx, omega_x[1,:] + omega_x[2,:] .== 0)
        @constraint(dcglp, cont[j=1:dim_t], omega_t[1,j] + omega_t[2,j] == 0) # must be in this form to recognize it as a vector

        assign_attributes!(dcglp, solver_param)
        
        disjunctiveCutsByIndex = [Vector{Hyperplane}() for i=1:data.dim_x]
        splits = Vector{Tuple{SparseVector{Float64, Int}, Float64}}()

        new(oracle_param, dcglp, typical_oracles, param, disjunctiveCutsByIndex, Vector{Hyperplane}(), splits)
    end
end