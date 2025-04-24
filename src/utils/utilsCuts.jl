export Hyperplane, aggregate
"""
hyperplane related structure and functions
"""
mutable struct Hyperplane
    
    a_x::SparseVector{Float64, Int} 
    a_t::SparseVector{Float64, Int} 
    a_0::Float64

    function Hyperplane(a_x::Vector{Float64}, 
        a_t::Vector{Float64},
        a_0::Float64)

        new(dropzeros!(sparsevec(a_x)), dropzeros!(sparsevec(a_t)), a_0)
    end

    function Hyperplane(dim_x::Int, dim_t::Int)
        # trivial hyperplane
        new(spzeros(dim_x), spzeros(dim_t), 0.0)
    end

    Hyperplane() = new()
end
"""
aggregate a set of hyperplanes into a single hyperplane via averaging
"""
function aggregate(hyperplanes::Vector{Hyperplane})
    h = Hyperplane()
    K = length(hyperplanes)
    h.a_x = sum(hyperplanes[j].a_x for j=1:K) * 1/K # averaged for numerical stability
    h.a_t = sum(hyperplanes[j].a_t for j=1:K) * 1/K
    h.a_0 = sum(hyperplanes[j].a_0 for j=1:K) * 1/K

    return h
end
"""
evaluate the violation of a hyperplane by (x_value, t_value)
"""
function evaluate_violation(h::Hyperplane, x_value::Vector{Float64}, t_value::Vector{Float64}; zero_tol = 1e-6)
    return h.a_0 + h.a_x' * x_value + h.a_t' * t_value >= zero_tol
end
"""
select a top fraction of a set of hyperplanes based on a score measured by a function f
"""
function select_top_fraction(a::Vector{Hyperplane}, f::Function, p::Float64)
    @assert 0 < p ≤ 1 "Fraction p must be in (0, 1]"
    
    # Apply function f to each element of a
    scores = f.(a)
    
    # Get the indices that would sort scores in descending order
    sorted_indices = sortperm(scores, rev=true)
    
    # How many elements to select
    k = ceil(Int, p * length(a))
    
    # Get the top-k indices and return corresponding elements from a
    top_indices = sorted_indices[1:k]
    return a[top_indices]
end
"""
parse a hyperplane to a JuMP expression
"""
function hyperplanes_to_expression(model::Model, hyperplanes::Vector{Hyperplane}, x_var::Vector{VariableRef}, t_var::Vector{VariableRef}, z_var::VariableRef)
    return @expression(model, [j in 1:length(hyperplanes)], hyperplanes[j].a_0 * z_var + hyperplanes[j].a_x' * x_var + hyperplanes[j].a_t' * t_var)
end
function hyperplanes_to_expression(model::Model, hyperplanes::Vector{Hyperplane}, x_var::Vector{VariableRef}, t_var::Vector{VariableRef})
    return @expression(model, [j in 1:length(hyperplanes)], hyperplanes[j].a_0 + hyperplanes[j].a_x' * x_var + hyperplanes[j].a_t' * t_var)
end

