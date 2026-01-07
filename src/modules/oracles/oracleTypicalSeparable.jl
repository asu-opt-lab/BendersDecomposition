export SeparableOracle, SeparableOracleParam

"""
Fallback constructor for subtypes of [`AbstractTypicalOracle`](@ref).

This method is invoked when the user attempts to construct a `SeparableOracle`
with a concrete oracle subtype `T` that does **not** define the required
type-call constructor:
```julia
    T(data::AbstractData, master::AbstractMaster;
      customize = customize_sub_model!, scen_idx::Int, param::AbstractOracleParam)
```
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

"""
    SeparableOracleParam <: AbstractOracleParam

Parameter type for [`SeparableOracle`](@ref).

`SeparableOracleParam` stores parameters specific to separable or
multi-scenario oracle aggregation. While it is currently empty, it exists to
support future extensions related to scenario handling, aggregation rules, or
termination logic.

See also: [`SeparableOracle`](@ref), [`AbstractOracleParam`](@ref)
"""
mutable struct SeparableOracleParam <: AbstractOracleParam
    # may contain parameters for scenario handling.
end

"""
    SeparableOracle <: AbstractTypicalOracle

Oracle wrapper for separable or multi-scenario Benders subproblems.

`SeparableOracle` constructs and manages a collection of independent
scenario-specific oracles, each respoSensible for generating cuts for a single
component of the auxiliary master variable `t`. 

This is particularly useful for:
- Multi-scenario stochastic programs
- Separable recourse functions

## Construction

```julia
SeparableOracle(
    data::AbstractData,
    master::Master,
    oracle::T,
    N::Int;
    customize = customize_sub_model!,
    sub_oracle_param::AbstractOracleParam = BasicOracleParam(),
    param::SeparableOracleParam = SeparableOracleParam(),
) where {T <: AbstractTypicalOracle}
```
#### Arguments
- `data`: User-defined problem data shared across all scenarios.
- `master`: The master module.
- `oracle`: A concrete oracle type `T` used to build each scenario-specific oracle.
- `N`: Number of scenarios (or separable components). This also defines the dimension of the master variable `t`.
- `customize`: User-provided function for building each subproblem model. It is called once per scenario with the corresponding `scen_idx`.
- `sub_oracle_param`: Parameters passed to each scenario-specific oracle.
- `param`: Parameters controlling the behavior of the `SeparableOracle` itself.

### Behavior
- Internally constructs `N` independent oracles of type `T`, each associated with a single scenario index `j = 1, …, N`.
- Each sub-oracle is assumed to generate cuts involving a scalar recourse variable `t[j]`.
- Cuts returned by sub-oracles are lifted to the full master space by embedding their `t`-coefficients in the appropriate component of the global `t` vector.

### Cut Generation
Calling `generate_cuts` on a `SeparableOracle`:
1. Delegates cut generation to each scenario-specific oracle.
2. Adjusts the dimensionality of the returned cuts to match the full master problem.
3. Aggregates all violated cuts across scenarios.
If all scenarios are satisfied, a trivial (zero) hyperplane is returned.

### Notes
- The concrete oracle type `T` must implement the constructor
    ```julia
    T(data::AbstractData, master::AbstractMaster;
    customize, scen_idx::Int, param::AbstractOracleParam)
    ```
    otherwise construction will fail with a descriptive error.
- Each scenario oracle is treated independently; no coupling across scenarios is assumed.
See also: `ClassicalOracle`, `generate_cuts`, `AbstractTypicalOracle`
"""
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






