# Base solve! function
function solve!(env::AbstractBendersEnv, solution_procedure::AbstractSolutionProcedure, cut_strategy::AbstractCutStrategy)
    throw(ArgumentError("Unsupported solution procedure type"))
end

# Solve with GenericCallback
function solve!(env::AbstractBendersEnv, solution_procedure::GenericCallback, cut_strategy::AbstractCutStrategy)
    function lazy_callback_wrapper(cb_data)
        solution_procedure.lazy_callback(cb_data, env, cut_strategy)
    end
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback_wrapper)
    
    if isnotnothing(solution_procedure.user_callback)
        function user_callback_wrapper(cb_data)
            solution_procedure.user_callback(cb_data, env, cut_strategy)
        end
        set_attribute(env.master.model, MOI.UserCutCallback(), user_callback_wrapper)
    end

    set_time_limit_sec(env.master.model, solution_procedure.time_limit)
    set_optimizer_attribute(env.master.model, MOI.Silent(), false)

    start_time = time()
    JuMP.optimize!(env.master.model)
    @info "node count" JuMP.node_count(env.master.model)
    @info "elapsed time" time() - start_time
    @info "objective bound" JuMP.objective_bound(env.master.model)
    @info "objective value" JuMP.objective_value(env.master.model)
    @info "relative gap" JuMP.relative_gap(env.master.model)
end

# Solve with CPLEX solver
function solve!(env::BendersEnv, solution_procedure::SolverCallback{CPLEX.Optimizer}, cut_strategy::CutStrategy)
    function callback_function_wrapper(cb_data::CPLEX.CallbackContext, context_id::Clong)
        solution_procedure.callback_function(cb_data, context_id, env, cut_strategy)
    end
    MOI.set(env.master.model, CPLEX.CallbackFunction(), callback_function_wrapper)
    set_time_limit_sec(env.master.model, solution_procedure.time_limit)
    set_optimizer_attribute(env.master.model, MOI.Silent(), false)
end

# Solve with Gurobi solver
function solve!(env::BendersEnv, solution_procedure::SolverCallback{Gurobi.Optimizer}, cut_strategy::CutStrategy)
    function callback_function_wrapper(cb_data, cb_where::Cint)
        solution_procedure.callback_function(cb_data, cb_where, env, cut_strategy)
    end
    MOI.set(env.master.model, Gurobi.CallbackFunction(), callback_function_wrapper)
    set_time_limit_sec(env.master.model, solution_procedure.time_limit)
    set_optimizer_attribute(env.master.model, MOI.Silent(), false)
end

# Base add_cuts! function
function add_cuts!(env::BendersEnv, expressions::Union{Vector{Any}, Any}, sub_obj_values::Union{Vector{Float64}, Float64}, cb_data)
    throw(ArgumentError("Unsupported types for add_cuts!"))
end

# Add cuts for vector inputs
function add_cuts!(env::BendersEnv, expressions::Vector{Any}, sub_obj_values::Vector{Float64}, cb_data)
    for (idx, (expr, sub_obj)) in enumerate(zip(expressions, sub_obj_values))
        if env.master.t_value[idx] <= sub_obj - 1e-06
            cut = @build_constraint(0 >= expr)
            MOI.submit(env.master.model, MOI.LazyConstraint(cb_data), cut)
        end
    end
end

# Add cuts for single expression
function add_cuts!(env::BendersEnv, expression::Any, sub_obj_value::Float64, cb_data)
    if env.master.t_value <= sub_obj_value - 1e-06
        cut = @build_constraint(0 >= expression)
        MOI.submit(env.master.model, MOI.LazyConstraint(cb_data), cut)
    end
end

# for base cut strategy
function root_node_preprocessing!(env::BendersEnv, solution_procedure::AbstractSequential, cut_strategy::AbstractCutStrategy)
    relax_integrality(env.master.model)
    df = solve!(env, solution_procedure, cut_strategy)
    set_binary.(env.master.model[:x])
    println("--------------------------------Root node preprocessing finished--------------------------------")
    return df
end

# function submit(model::ModelLike, sub::AbstractSubmittable, args...)
#     if supports(model, sub)
#         throw(
#             ArgumentError(
#                 "Submitting $(typeof.(args)) for `$(typeof(sub))` is not valid.",
#             ),
#         )
#     else
#         throw(
#             UnsupportedSubmittable(
#                 sub,
#                 "submit(::$(typeof(model)), ::$(typeof(sub))) is not supported.",
#             ),
#         )
#     end
# end