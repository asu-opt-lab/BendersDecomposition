export SeparableOracle, SeparableOracleParam

mutable struct SeparableOracleParam <: AbstractOracleParam
    # may contain parameters for scenario handling.
end

mutable struct SeparableOracle <: AbstractTypicalOracle
    oracle_param::SeparableOracleParam 

    oracles::Vector{AbstractTypicalOracle}
    N::Int

    function SeparableOracle(data::Data, 
                            oracle::T, 
                            N::Int; 
                            solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9), 
                            sub_oracle_param::AbstractOracleParam = EmptyOracleParam(),
                            oracle_param::SeparableOracleParam = SeparableOracleParam()) where {T<:AbstractTypicalOracle}
        @debug "Building classical separable oracle"
        # assume each oracle is associated with a single t, that is dim_t = N
        oracles = typeof(sub_oracle_param) != EmptyOracleParam ? [T(data, scen_idx=j, solver_param = solver_param, oracle_param = sub_oracle_param) for j=1:N] : [T(data, scen_idx=j, solver_param = solver_param) for j=1:N]

        new(oracle_param, oracles, N)
    end
end

function generate_cuts(oracle::SeparableOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol = 1e-6, time_limit = 3600)
    tic = time()
    N = oracle.N
    is_in_L = Vector{Bool}(undef,N)
    sub_obj_val = Vector{Vector{Float64}}(undef,N)
    hyperplanes = Vector{Vector{Hyperplane}}(undef,N)
    
    # do threads?
    for j=1:N
        is_in_L[j], hyperplanes[j], sub_obj_val[j] = generate_cuts(oracle.oracles[j], x_value, [t_value[j]]; tol=tol, time_limit=get_sec_remaining(tic, time_limit))

        # correct dimension for t_j's
        for h in hyperplanes[j]
            coeff_t = h.a_t[1]
            h.a_t = spzeros(length(t_value)) 
            h.a_t[j] = coeff_t
        end
    end

    if any(.!is_in_L)
        return false, reduce(vcat, hyperplanes), reduce(vcat, sub_obj_val)
    else
        return true, [Hyperplane(length(x_value), length(t_value))], t_value
    end
end






