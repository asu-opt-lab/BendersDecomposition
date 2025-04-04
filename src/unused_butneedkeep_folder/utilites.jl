"""
Solve the master problem
"""
function solve_master!(master::AbstractMasterProblem)
    optimize!(master.model)
    if termination_status(master.model) == MOI.OPTIMAL
        master.objective_value = JuMP.objective_value(master.model)
        master.integer_variable_values = value.(master.var[:x])
        master.continuous_variable_values = value.(master.var[:t])
    else
        throw(ErrorException("Master problem not solved to optimality"))
    end
end

"""
Solve the sub problem
"""
function solve_sub!(sub::AbstractSubProblem, x_value::Vector{Float64})
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

"""
Check termination criteria based on gap, time limit and iteration limit
"""
# 所有参数都是 nothing
function is_terminated(state::BendersState, ::Nothing, ::Nothing, ::Nothing)
    return state.gap <= 0.001
end

# 只有 iteration_limit 不是 nothing
function is_terminated(state::BendersState, iteration_limit::Int, ::Nothing, ::Nothing)
    return state.iteration >= iteration_limit || state.gap <= 0.001
end

# 只有 time_limit 不是 nothing
function is_terminated(state::BendersState, ::Nothing, time_limit::Float64, ::Nothing)
    return state.total_time >= time_limit || state.gap <= 0.001
end

# 只有 gap_tolerance 不是 nothing
function is_terminated(state::BendersState, ::Nothing, ::Nothing, gap_tolerance::Float64)
    return state.gap <= gap_tolerance || state.gap <= 0.001
end

# iteration_limit 和 time_limit 不是 nothing
function is_terminated(state::BendersState, iteration_limit::Int, time_limit::Float64, ::Nothing)
    return state.iteration >= iteration_limit || state.total_time >= time_limit
end

# iteration_limit 和 gap_tolerance 不是 nothing
function is_terminated(state::BendersState, iteration_limit::Int, ::Nothing, gap_tolerance::Float64)
    return state.iteration >= iteration_limit || state.gap <= gap_tolerance 
end

# time_limit 和 gap_tolerance 不是 nothing
function is_terminated(state::BendersState, ::Nothing, time_limit::Float64, gap_tolerance::Float64)
    return state.total_time >= time_limit || state.gap <= gap_tolerance
end

# 所有参数都不是 nothing
function is_terminated(state::BendersState, iteration_limit::Int, time_limit::Float64, gap_tolerance::Float64)
    return state.iteration >= iteration_limit || 
           state.total_time >= time_limit || 
           state.gap <= gap_tolerance
end

mutable struct BendersState
    iter::Int
    LB::Float64
    UB::Float64
    gap::Float64
    master_time::Float64
    sub_time::Float64
    total_time::Float64

    function BendersState()
        return new(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
end

function add_cuts!(env::BendersEnv, cut::JuMP.ConstraintRef)
    @constraint(env.master.model, 0 >= cut)
end

function add_cuts!(env::BendersEnv, cuts::Vector{JuMP.ConstraintRef})
    for cut in cuts
        add_cut!(env, cut)
    end
end
