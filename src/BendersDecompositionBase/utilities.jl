mutable struct BendersState
    iteration::Int
    LB::Float64
    UB::Float64
    gap::Float64
    master_time::Float64
    sub_time::Float64
    total_time::Float64

    function BendersState()
        return new(0, -Inf, Inf, Inf, 0.0, 0.0, 0.0)
    end
end


"""
Solve the master problem
"""
function solve_master!(master::AbstractMasterProblem)
    optimize!(master.model)
    if termination_status(master.model) == MOI.OPTIMAL
        master.objective_value = JuMP.objective_value(master.model)
        master.integer_variable_values = value.(master.variables[:integer_variables])
        master.continuous_variable_values = value.(master.variables[:continuous_variables])
    else
        throw(ErrorException("Master problem not solved to optimality, status: $(termination_status(master.model))"))
    end
end

"""
Solve the sub problem
"""
function solve_sub!(sub::AbstractSubProblem, x_value::Vector{Float64})
    sub.fixed_x_values = x_value
    set_normalized_rhs.(sub.fixed_x_constraints, x_value)
    optimize!(sub.model)
end

"""
Print iteration information if verbose mode is on
"""
function print_iteration_info(state::BendersState)
    @printf("Iter: %4d | LB: %12.4f | UB: %11.4f | Gap: %8.3f%% | Time: (M: %6.2f, S: %6.2f) | Elapsed: %6.2f\n",
           state.iteration, state.LB, state.UB, state.gap, 
           state.master_time, state.sub_time, state.total_time)
end

function add_cuts!(env::AbstractBendersEnv, cut::AffExpr)
    @constraint(env.master.model, 0 >= cut)
    @info "Added cut: $cut"
end

function add_cuts!(env::AbstractBendersEnv, cuts::Vector{AffExpr})
    for cut in cuts
        add_cuts!(env, cut)
    end
end


function update_upper_bound_and_gap!(state::BendersState, env::AbstractBendersEnv, sub_obj_val::Float64)

    state.UB = min(state.UB, sub_obj_val + env.master.objective_value - env.master.continuous_variable_values)

    update_gap!(state)
end

function update_gap!(state::BendersState)
    state.gap = (state.UB - state.LB) / state.UB * 100
end

"""
    is_terminated(state::BendersState, iteration_limit::Union{Int, Nothing}=nothing, 
                 time_limit::Union{Float64, Nothing}=nothing, 
                 gap_tolerance::Union{Float64, Nothing}=0.001)

Check if the Benders decomposition algorithm should terminate based on the given criteria.
Returns true if any of the following conditions are met:
- The iteration count exceeds iteration_limit (if provided)
- The total time exceeds time_limit (if provided)
- The gap is less than or equal to gap_tolerance (defaults to 0.001 if not provided)
"""
function is_terminated(
    state::BendersState,
    iteration_limit::Union{Int, Nothing}=nothing,
    time_limit::Union{Float64, Nothing}=nothing,
    gap_tolerance::Union{Float64, Nothing}=0.01
)
    if !isnothing(iteration_limit) && state.iteration >= iteration_limit
        return true
    end
    
    if !isnothing(time_limit) && state.total_time >= time_limit
        return true
    end
    
    if !isnothing(gap_tolerance) && state.gap <= gap_tolerance
        return true
    end

    return false
end