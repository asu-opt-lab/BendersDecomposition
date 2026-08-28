"""
    (::Type{T})(data::AbstractData, master::AbstractMaster;
                model = update_sub_model!,
                scen_idx::Int,
                param::AbstractOracleParam,
                optimizer = DEFAULT_OPTIMIZER) where T <: AbstractTypicalOracle

Fallback constructor for [`AbstractTypicalOracle`](@ref) subtypes used by [`SeparableOracle`](@ref).

A concrete typical-oracle type intended for use with `SeparableOracle` must provide a constructor with this interface so that `SeparableOracle` can instantiate one sub-oracle for each scenario.

The `model` keyword specifies the subproblem model-update function, and `scen_idx` identifies the scenario or subproblem associated with the oracle.

# Throws

Throws an error if `T` does not implement the required constructor.
"""
(::Type{T})(data::AbstractData, master::AbstractMaster;
            model = update_sub_model!,
            scen_idx::Int,
            param::AbstractOracleParam,
            optimizer = DEFAULT_OPTIMIZER) where T <: AbstractTypicalOracle =
    throw(UnimplementedInterfaceException(
        """
        SeparableOracle: 
        Oracle subtype $(T) does not implement the required constructor needed by `SeparableOracle`.

        Expected constructor signature:

          $(T)(data::AbstractData, master::AbstractMaster;
              model = update_sub_model!, scen_idx::Int, param::AbstractOracleParam,
              optimizer = ...)

        Define this constructor for $(T) in order to use it with `SeparableOracle`.
        """
    ))

"""
    SeparableOracleParam <: AbstractOracleParam

Parameters controlling [`SeparableOracle`](@ref).

This parameter container is currently empty and serves as an extension point for future controls related to scenario handling or parallel evaluation.

See also: [`SeparableOracle`](@ref), [`AbstractOracleParam`](@ref)
"""
mutable struct SeparableOracleParam <: AbstractOracleParam
    # may contain parameters for scenario handling.
end

"""
    SeparableOracle <: AbstractTypicalOracle

Composite oracle for problems with multiple independent subproblems.

`SeparableOracle` constructs one [`AbstractTypicalOracle`](@ref) of the specified concrete type for each scenario or subproblem and evaluates the sub-oracles in parallel. Each sub-oracle generates cuts associated with one component of the auxiliary variable `t`.

The oracle type is specified by `oracle`, whose concrete type must provide the constructor required to instantiate a sub-oracle for each scenario.

# Fields
param::SeparableOracleParam: Parameters controlling the separable oracle.
oracles::Vector{AbstractTypicalOracle}: Typical sub-oracles, one for each scenario or subproblem.
N::Int: Number of subproblems.

# Constructor

    SeparableOracle(
        data::AbstractData,
        master::Master,
        oracle_type::Type{T},
        N::Int;
        model = update_sub_model!,
        sub_oracle_param::AbstractOracleParam = BasicOracleParam(),
        param::SeparableOracleParam = SeparableOracleParam(),
        optimizer = DEFAULT_OPTIMIZER,
    ) where T <: AbstractTypicalOracle


Construct a separable oracle with `N` sub-oracles of type `oracle_type`.

A new sub-oracle is constructed for each scenario using `scen_idx = 1:N`. The same `sub_oracle_param` and `model` function are passed to every sub-oracle.

`N` must equal `master.dim_t`, with one component of the master auxiliary variable `t` associated with each subproblem.

# Throws

Throws a `DimensionMismatch` if `N != master.dim_t`.

See also: [`ClassicalOracle`](@ref), [`ParetoOracle`](@ref), [`UnifiedOracle`](@ref)
"""
mutable struct SeparableOracle <: AbstractTypicalOracle
    param::SeparableOracleParam 

    oracles::Vector{AbstractTypicalOracle}
    N::Int

    function SeparableOracle(data::AbstractData, 
                            master::Master,
                            oracle_type::Type{T}, 
                            N::Int; 
                            model = update_sub_model!,
                            sub_oracle_param::AbstractOracleParam = BasicOracleParam(),
                            param::SeparableOracleParam = SeparableOracleParam(),
                            optimizer = DEFAULT_OPTIMIZER) where {T<:AbstractTypicalOracle}
        @debug "Building classical separable oracle"
        @info "SeparableOracle: N=$N subproblems, $(Threads.nthreads()) threads available for parallel execution"
        
        N == master.dim_t || throw(DimensionMismatch("SeparableOracle: `N` must equal master.dim_t ($(master.dim_t)), got $N."))

        oracles = [oracle_type(data, master; model = model, scen_idx = j, param = sub_oracle_param, optimizer = optimizer) for j in 1:N]

        new(param, oracles, N)
    end
end

"""
    generate_cuts(
        oracle::SeparableOracle,
        x_value::Vector{Float64},
        t_value::Vector{Float64};
        tol_normalize = 1.0,
        time_limit = 3600.0,
    )

Generate Benders cuts by evaluating all sub-oracles in parallel.

Each sub-oracle is evaluated at the common candidate `x_value` and its corresponding component `t_value[j]`. Generated cuts are expanded to the full `t` dimension and associated with the corresponding subproblem.

If any sub-oracle separates the candidate, the generated cuts and subproblem objective values are returned collectively. Otherwise, the candidate is reported as belonging to the separable oracle's feasible region.
"""
function generate_cuts(oracle::SeparableOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600.0)
    tic = time()
    N = oracle.N
    is_in_L = Vector{Bool}(undef,N)
    sub_obj_val = Vector{Vector{Float64}}(undef,N)
    hyperplanes = Vector{Vector{Hyperplane}}(undef,N)

    try
        Threads.@threads for j=1:N
            is_in_L[j], hyperplanes[j], sub_obj_val[j] = generate_cuts(oracle.oracles[j], x_value, [t_value[j]], tol_normalize = tol_normalize; time_limit = get_sec_remaining(tic, time_limit))

            # correct dimension for t_j's
            for h in hyperplanes[j]
                coeff_t = h.a_t[1]
                h.a_t = spzeros(length(t_value))
                h.a_t[j] = coeff_t
            end
        end
    catch err
        err isa CompositeException || rethrow()
        task_failure = first(err.exceptions)
        task_failure isa TaskFailedException || rethrow()
        failure = first(current_exceptions(task_failure.task; backtrace = false))
        throw(failure.exception)
    end

    if any(.!is_in_L)
        return false, reduce(vcat, hyperplanes), reduce(vcat, sub_obj_val)
    else
        return true, [Hyperplane(length(x_value), length(t_value))], deepcopy(t_value)
    end
end
