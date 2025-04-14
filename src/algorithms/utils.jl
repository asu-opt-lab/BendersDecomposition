export get_sec_remaining

"""
Print iteration information if verbose mode is on
"""
function print_iteration_info(state::AbstractLoopState, log::AbstractLoopLog)
    throw(UndefError("ddd"))
end

function update_upper_bound_and_gap!(state::AbstractLoopState, f::Function)
    evaluation = f(state.f_x, state.x_value)
    state.UB = min(state.UB, evaluation)
    state.gap = (state.UB - state.LB) / abs(state.UB) * 100
end
function is_terminated(state::AbstractLoopState, log::AbstractLoopLog, params::BendersParams)
    throw(UndefError("ddd33"))
end

function check_lb_improvement!(state::AbstractLoopState, log::AbstractLoopLog; zero_tol = 1e-8, tol_imprv = 1e-4)
    prev_lb = log.iterations[end].LB
    lb_improvement = abs(prev_lb) < zero_tol ? abs(state.LB - prev_lb) : abs((state.LB - prev_lb) / prev_lb) * 100
    
    # Check for improvement
    if lb_improvement < tol_imprv
        log.consecutive_no_improvement += 1
    else
        # Reset counter if there's improvement
        log.consecutive_no_improvement = 0
    end
end

function add_constraints(model::Model, constr_symbol::Symbol, exprs::Vector{AffExpr})
    # add constraints in the form of 0 .>= expr
    if haskey(model, constr_symbol)
        append!(model[constr_symbol], @constraint(model, 0 .>= exprs))
    else
        model[constr_symbol] = @constraint(model, 0 .>= exprs)
    end
end

function evaluate_violation(h::Hyperplane, x_value::Vector{Float64}, t_value::Vector{Float64}; zero_tol = 1e-6)
    return h.a_0 + h.a_x' * x_value + h.a_t' * t_value >= zero_tol
end

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

function hyperplanes_to_expression(model::Model, hyperplanes::Vector{Hyperplane}, x_var::Vector{VariableRef}, t_var::Vector{VariableRef}, z_var::VariableRef)
    return @expression(model, [j in 1:length(hyperplanes)], hyperplanes[j].a_0 * z_var + hyperplanes[j].a_x' * x_var + hyperplanes[j].a_t' * t_var)
end
function hyperplanes_to_expression(model::Model, hyperplanes::Vector{Hyperplane}, x_var::Vector{VariableRef}, t_var::Vector{VariableRef})
    return @expression(model, [j in 1:length(hyperplanes)], hyperplanes[j].a_0 + hyperplanes[j].a_x' * x_var + hyperplanes[j].a_t' * t_var)
end

# Helper functions for Benders DecompositionLog
function get_sec_remaining(tic::Float64, time_limit::Float64; tol = 1e-4)
    # tol needed to prevent parameter being too small
    time_elapsed = time() - tic
    return max(time_limit - time_elapsed, tol)
end
function get_sec_remaining(log::AbstractLoopLog, param::BendersParams)
    return get_sec_remaining(log.start_time, param.time_limit)
end
function record_iteration!(log::AbstractLoopLog, state::AbstractLoopState)
    push!(log.iterations, state)
    log.n_iter += 1
end

function to_dataframe(log::AbstractLoopLog)
    return DataFrame(
        LB = [info.LB for info in log.iterations],
        UB = [info.UB for info in log.iterations],
        gap = [info.gap for info in log.iterations],
        master_time = [info.master_time for info in log.iterations],
        oracle_time = [info.oracle_time for info in log.iterations],
        total_time = [info.total_time for info in log.iterations],
        # is_in_L = [info.is_in_L for info in log.iterations]
    )
end

