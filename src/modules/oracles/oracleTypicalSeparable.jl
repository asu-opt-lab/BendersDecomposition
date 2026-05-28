
"""
Fallback constructor for subtypes of [`AbstractTypicalOracle`](@ref).

This method is invoked when the user attempts to construct a `SeparableOracle`
with a concrete oracle subtype `T` that does **not** define the required
type-call constructor:

    T(data::AbstractData, master::AbstractMaster;
      model = customize_sub_model!, scen_idx::Int, param::AbstractOracleParam)

Calling this fallback indicates that the oracle type `T` has not implemented
the interface expected by `SeparableOracle`. Any concrete oracle intended for
use with `SeparableOracle` must therefore define a constructor matching the
signature above.

# Errors
Throws an error indicating that the subtype `T` must provide the required
constructor.
"""
(::Type{T})(data::AbstractData, master::AbstractMaster;
            model = customize_sub_model!,
            customize = nothing,
            scen_idx::Int,
            param::AbstractOracleParam,
            optimizer = DEFAULT_OPTIMIZER) where T <: AbstractTypicalOracle =
    throw(UndefError(
        """
        Oracle subtype $(T) does not implement the required constructor
        needed by `SeparableOracle`.

        Expected constructor signature:

          $(T)(data::AbstractData, master::AbstractMaster;
              model = customize_sub_model!, scen_idx::Int,
              param::AbstractOracleParam, optimizer = ...)

        Define this constructor for $(T) in order to use it with `SeparableOracle`.
        """
    ))

"""
    SeparableOracleParam <: AbstractOracleParam

Parameter container for [`SeparableOracle`](@ref).

The type is currently empty and acts as a placeholder for future controls over
scenario aggregation, scheduling, or parallel evaluation.

See also: [`SeparableOracle`](@ref), [`AbstractOracleParam`](@ref)
"""
mutable struct SeparableOracleParam <: AbstractOracleParam
    # may contain parameters for scenario handling.
end

"""
    SeparableOracle <: AbstractTypicalOracle

Wrapper oracle for problems with multiple independent recourse subproblems.

`SeparableOracle` builds one typical oracle per scenario or block and evaluates
them in parallel. Each sub-oracle is responsible for one component of the
vector-valued recourse approximation `t`.

# Constructor
```julia
SeparableOracle(data, master, oracle_template, N;
                model = customize_sub_model!,
                sub_oracle_param = BasicOracleParam(),
                param = SeparableOracleParam())
```

Here `oracle_template` is any instance whose type is a subtype of
`AbstractTypicalOracle`; `SeparableOracle` uses that type to instantiate one
sub-oracle per scenario. `N` is the number of subproblems.

See also: [`ClassicalOracle`](@ref), [`ParetoOracle`](@ref), [`UnifiedOracle`](@ref)
"""
mutable struct SeparableOracle <: AbstractTypicalOracle
    param::SeparableOracleParam 

    oracles::Vector{AbstractTypicalOracle}
    N::Int

    function SeparableOracle(data::AbstractData, 
                            master::Master,
                            oracle::T, 
                            N::Int; 
                            model = customize_sub_model!,
                            customize = nothing,
                            sub_oracle_param::AbstractOracleParam = BasicOracleParam(),
                            param::SeparableOracleParam = SeparableOracleParam(),
                            optimizer = DEFAULT_OPTIMIZER) where {T<:AbstractTypicalOracle}
        @debug "Building classical separable oracle"
        @info "SeparableOracle: N=$N subproblems, $(Threads.nthreads()) threads available for parallel execution"
        model_update = _resolve_model_update_keyword(model, customize)
        # assume each oracle is associated with a single t, that is dim_t = N
        oracles = [T(data, master; model = model_update, scen_idx = j, param = sub_oracle_param, optimizer = optimizer) for j in 1:N]

        new(param, oracles, N)
    end
end

function generate_cuts(oracle::SeparableOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600.0)
    tic = time()
    N = oracle.N
    is_in_L = Vector{Bool}(undef,N)
    sub_obj_val = Vector{Vector{Float64}}(undef,N)
    hyperplanes = Vector{Vector{Hyperplane}}(undef,N)

    Threads.@threads for j=1:N
        is_in_L[j], hyperplanes[j], sub_obj_val[j] = generate_cuts(oracle.oracles[j], x_value, [t_value[j]], tol_normalize = tol_normalize; time_limit = get_sec_remaining(tic, time_limit))

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
        return true, [Hyperplane(length(x_value), length(t_value))], deepcopy(t_value)
    end
end
