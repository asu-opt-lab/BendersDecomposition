export get_sec_remaining

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
function get_sec_remaining(log::BendersLog, param::BendersParams)
    return get_sec_remaining(log.start_time, param.time_limit)
end
function record_iteration!(log::BendersLog, state::BendersState)
    push!(log.iterations, state)
    log.LB = state.LB
    log.UB = state.UB
end

function update_upper_bound_and_gap!(state::BendersState, env::BendersEnv, sub_obj_val::Vector{Float64})
    state.UB = min(state.UB, env.data.c_t' * sub_obj_val + env.data.c_x' * env.master.x_value)
    state.gap = (state.UB - state.LB) / abs(state.UB) * 100
end

function to_dataframe(log::BendersLog)
    return DataFrame(
        iter = [info.iteration for info in log.iterations],
        LB = [info.LB for info in log.iterations],
        UB = [info.UB for info in log.iterations],
        gap = [info.gap for info in log.iterations],
        master_time = [info.master_time for info in log.iterations],
        oracle_time = [info.oracle_time for info in log.iterations],
        total_time = [info.total_time for info in log.iterations],
        is_in_L = [info.is_in_L for info in log.iterations]
    )
end

