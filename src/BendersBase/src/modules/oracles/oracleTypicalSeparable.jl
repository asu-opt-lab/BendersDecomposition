export SeparableOracle, SeparableOracleParam

"""
Fallback constructor for subtypes of [`AbstractTypicalOracle`](@ref).

This method is invoked when the user attempts to construct a `SeparableOracle`
with a concrete oracle subtype `T` that does **not** define the required
type-call constructor:

    T(data::AbstractData, master::AbstractMaster;
      customize = customize_sub_model!, scen_idx::Int, param::AbstractOracleParam)

Calling this fallback indicates that the oracle type `T` has not implemented
the interface expected by `SeparableOracle`. Any concrete oracle intended for
use with `SeparableOracle` must therefore define a constructor matching the
signature above.

# Errors
Throws an error indicating that the subtype `T` must provide the required
constructor.
"""
(::Type{T})(data::AbstractData, master::AbstractMaster;
            customize = customize_sub_model!,
            scen_idx::Int,
            param::AbstractOracleParam) where T <: AbstractTypicalOracle =
    throw(UndefError(
        """
        Oracle subtype $(T) does not implement the required constructor
        needed by `SeparableOracle`.

        Expected constructor signature:

          $(T)(data::AbstractData, master::AbstractMaster;
              customize = customize_sub_model!, scen_idx::Int, param::AbstractOracleParam)

        Define this constructor for $(T) in order to use it with `SeparableOracle`.
        """
    ))


mutable struct SeparableOracleParam <: AbstractOracleParam
    # may contain parameters for scenario handling.
end

mutable struct SeparableOracle <: AbstractTypicalOracle
    param::SeparableOracleParam 

    oracles::Vector{AbstractTypicalOracle}
    N::Int

    function SeparableOracle(data::AbstractData, 
                            master::Master,
                            oracle::T, 
                            N::Int; 
                            customize = customize_sub_model!,
                            sub_oracle_param::AbstractOracleParam = BasicOracleParam(),
                            param::SeparableOracleParam = SeparableOracleParam()) where {T<:AbstractTypicalOracle}
        @debug "Building classical separable oracle"
        # assume each oracle is associated with a single t, that is dim_t = N
        oracles = [T(data, master; customize = customize, scen_idx=j, param = sub_oracle_param) for j=1:N] 

        new(param, oracles, N)
    end
end

function generate_cuts(oracle::SeparableOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600.0)
    tic = time()
    N = oracle.N
    is_in_L = Vector{Bool}(undef,N)
    sub_obj_val = Vector{Vector{Float64}}(undef,N)
    hyperplanes = Vector{Vector{Hyperplane}}(undef,N)

    for j=1:N
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






