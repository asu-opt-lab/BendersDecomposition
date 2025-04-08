export ClassicalOracle, SeparableOracle, DisjunctiveOracle, generate_cut_coefficients, Hyperplane, aggregate


mutable struct Hyperplane
    
    a_x::Vector{Float64}
    a_t::Vector{Float64}
    a_0::Float64

    # a_x_indices
    # a_x_value
    function Hyperplane(a_x::Vector{Float64}, 
        a_t::Vector{Float64},
        a_0::Float64)
        
        new(a_x, a_t, a_0)
    end

    function Hyperplane(dim_x::Int, dim_t::Int)
        # trivial hyperplane
        new(zeros(dim_x), zeros(dim_t), 0.0)
    end

    Hyperplane() = new()
end

function aggregate(hyperplanes::Vector{Hyperplane})
    h = Hyperplane()
    K = length(hyperplanes)
    h.a_x = sum(hyperplanes[j].a_x for j=1:K)
    h.a_t = sum(hyperplanes[j].a_t for j=1:K)
    h.a_0 = sum(hyperplanes[j].a_0 for j=1:K)

    return h
end

mutable struct ClassicalOracle <: AbstractTypicalOracle
    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    other_constraints::Vector{ConstraintRef}

    function ClassicalOracle(data::Data; scen_idx::Int=-1)
        @debug "Building classical oracle"
        model = Model()

        # Define coupling variables and constraints
        @variable(model, x[1:data.dim_x])
        @constraint(model, fix_x, x .== 0)

        other_constr = Vector{ConstraintRef}()

        new(model, fix_x, other_constr)
    end

    ClassicalOracle() = new()
end

mutable struct SeparableOracle <: AbstractTypicalOracle
    oracles::Vector{AbstractTypicalOracle}
    N::Int

    function SeparableOracle(data::Data, oracle::T, N::Int) where {T<:AbstractTypicalOracle}
        @debug "Building classical separable oracle"
        # assume each oracle is associated with a single t, that is dim_t = N
        oracles = [T(data, scen_idx=j) for j=1:N]

        new(oracles, N)
    end
end

function generate_cuts(oracle::ClassicalOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol = 1e-6)
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    optimize!(oracle.model)
    status = dual_status(oracle.model)

    h = Hyperplane(length(x_value), length(t_value))

    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(oracle.model)

        if sub_obj_val >= t_value[1] + tol
            h.a_x = dual.(oracle.fixed_x_constraints) 
            h.a_t = [-1.0] 
            h.a_0 = sub_obj_val - h.a_x'*x_value 
            
            return false, [h], [sub_obj_val]
        end
        
        return true, [h], t_value

    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(oracle.model)
            h.a_x = dual.(oracle.fixed_x_constraints)
            h.a_t = [0.0]
            h.a_0 = dual.(oracle.other_constraints)'*normalized_rhs.(oracle.other_constraints)
            # @info coefficients_t' * t_value + coefficients_x' * x_value + constant_term
            return false, [h], [Inf]
        else
            throw(ErrorException("ClassicalOracle: Infeasible subproblem has no dual solution"))
        end
        
    else
        throw(ErrorException("ClassicalOracle: Unexpected dual status"))
    end
end

function generate_cuts(oracle::SeparableOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol = 1e-6)
    N = oracle.N
    is_in_L = Vector{Bool}(undef,N)
    sub_obj_val = Vector{Vector{Float64}}(undef,N)
    hyperplanes = Vector{Vector{Hyperplane}}(undef,N)
    
    for j=1:N
        is_in_L[j], hyperplanes[j], sub_obj_val[j] = generate_cuts(oracle.oracles[j], x_value, [t_value[j]]; tol=tol)

        # correct dimension for t_j's
        for h in hyperplanes[j]
            coeff_t = h.a_t[1]
            h.a_t = zeros(length(t_value)) 
            h.a_t[j] = coeff_t
        end
    end

    if any(.!is_in_L)
        return false, reduce(vcat, hyperplanes), reduce(vcat, sub_obj_val)
    else
        return true, [Hyperplane(length(x_value), length(t_value))], t_value
    end
end






