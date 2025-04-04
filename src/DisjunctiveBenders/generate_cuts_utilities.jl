mutable struct DCGLPState
    iteration::Int
    LB::Float64
    UB::Float64
    UB_k::Float64
    UB_v::Float64
    gap::Float64
    master_time::Float64
    sub_k_time::Float64
    sub_v_time::Float64
    total_time::Float64

    function DCGLPState()
        return new(0, -Inf, Inf, Inf, Inf, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
end


function solve_and_get_dcglp_values(model::Model, ::LNorm)
    optimize!(model)
    k_values = (constant = value(model[:k₀]), x = value.(model[:kₓ]), t = value.(model[:kₜ]))
    v_values = (constant = value(model[:v₀]), x = value.(model[:vₓ]), t = value.(model[:vₜ]))
    other_values = (τ = value(model[:τ]), sₓ = value.(model[:sₓ]), sₜ = value.(model[:sₜ]))
    return k_values, v_values, other_values
end

function solve_and_get_cut_coefficients(sub::AbstractSubProblem, korv_values, base_cut_strategy::AbstractCutStrategy)
    if isapprox(korv_values.constant, 0.0, atol=1e-05)
        return  [], max.(0, korv_values.t)
    end
    input_x = @. abs(korv_values.x / korv_values.constant) # abs value
    solve_sub!(sub, input_x)
    dual_info, obj_value = generate_cut_coefficients(sub, base_cut_strategy)
    obj_value *= korv_values.constant
    dual_info, obj_value = correct_cut_and_obj_values(dual_info, obj_value, korv_values.t) # remove redundant cuts
    return  dual_info, obj_value
end

"""
    correct_cut_and_obj_values(dual_info::Vector, obj_value::Float64, t_values::Float64)

Correct the cut coefficients and objective value for a single scenario.
Returns empty dual_info and t_values if the objective value is close to t_values.
"""
function correct_cut_and_obj_values(dual_info::Any, obj_value::Float64, t_values::Float64)
    if obj_value <= t_values + 1e-04
        return [], t_values
    end
    return dual_info, obj_value
end

"""
    correct_cut_and_obj_values(dual_info::Vector, obj_value::Vector{Float64}, t_values::Vector{Float64})

Correct the cut coefficients and objective values for multiple scenarios.
Removes cuts where the objective value is close to t_values.
"""
function correct_cut_and_obj_values(dual_info::Any, obj_value::Vector{Float64}, t_values::Vector{Float64})
    valid_indices = obj_value .> t_values .+ 1e-04
    obj_value[.!valid_indices] .= t_values[.!valid_indices]
    
    if !any(valid_indices)
        return [], obj_value
    end
    return dual_info[valid_indices], obj_value
end


function update_bounds!(state, k_values, v_values, other_values, obj_value_k, obj_value_v, t_value, norm_type::LNorm)
    state.LB = other_values.τ
    diff_st = obj_value_k .+ obj_value_v .- t_value
    state.UB = update_UB!(state.UB, other_values.sₓ, diff_st, norm_type)
    state.UB_k = obj_value_k .- k_values.t
    state.UB_v = obj_value_v .- v_values.t
    state.gap = (state.UB - state.LB)/abs(state.UB) * 100
end

update_UB!(UB,_sx,diff_st,::L1Norm) = min(UB,norm([ _sx; diff_st], Inf))
update_UB!(UB,_sx,diff_st,::L2Norm) = min(UB,norm([ _sx; diff_st], 2))
update_UB!(UB,_sx,diff_st,::LInfNorm) = min(UB,norm([ _sx; diff_st], 1))

function print_dcglp_iteration_info(state)
    @printf("   Iter: %4d | LB: %8.4f | UB: %8.4f | Gap: %6.2f%% | UB_k: %8.2f | UB_v: %8.2f | Master time: %6.2f | Sub_k time: %6.2f | Sub_v time: %6.2f \n",
           state.iteration, state.LB, state.UB, state.gap, sum(state.UB_k), sum(state.UB_v), state.master_time, state.sub_k_time, state.sub_v_time)
end

function is_terminated(
    state::DCGLPState,
    iteration_limit::Union{Int, Nothing}=30,
    time_limit::Union{Float64, Nothing}=200,
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
    
    if state.UB - state.LB <= 1e-03
        return true
    end
    
    return false
end

"""
    add_cuts_to_model!(env::BendersEnv, cuts::Union{AffExpr, Vector{AffExpr}}, constraints::Vector{ConstraintRef})

Add cuts to the model and store them in the constraints vector.
If cuts is a vector, each cut will be added as a separate constraint.
"""
function add_cuts_to_model!(env::DisjunctiveBendersEnv, cuts::Union{AffExpr, Vector{AffExpr}}, constraints::Vector{ConstraintRef})
    if cuts isa Vector
        for cut in cuts
            push!(constraints, @constraint(env.dcglp.model, 0 >= cut))
        end
    else
        push!(constraints, @constraint(env.dcglp.model, 0 >= cuts))
    end
end

"""
    add_cuts_with_strategy!(env::BendersEnv, dual_info::Tuple, cut_strategy::CutStrategy, primary_cuts, secondary_cuts)

Add cuts to the model based on the cut strategy.
"""
function add_cuts_with_strategy!(
    env::DisjunctiveBendersEnv,
    dual_info::Tuple,
    cut_strategy::DisjunctiveCut,
    primary_cuts,
    secondary_cuts
)
    cuts_k, cuts_v, cuts_master = build_dcglp_cuts(env, dual_info, cut_strategy.base_cut_strategy)
    
    # Add primary cuts
    add_cuts_to_model!(env, primary_cuts(cuts_k, cuts_v), env.dcglp.dcglp_constraints)
    
    # Add two-sided cuts if enabled
    if cut_strategy.use_two_sided_cuts 
        add_cuts_to_model!(env, secondary_cuts(cuts_k, cuts_v), env.dcglp.dcglp_constraints)
    end
    
    # Add master cuts if enabled
    if cut_strategy.include_master_cuts 
        append!(env.dcglp.master_cuts, cuts_master)
    end
end

"""
    CutType

Enum type to specify which type of cuts to add (k or v).
"""
@enum CutType begin
    K_CUTS
    V_CUTS
end

"""
    add_cuts!(env::BendersEnv, dual_info::Tuple, cut_strategy::CutStrategy, cut_type::CutType)

Add cuts to the model based on the cut type (k or v) and strategy.
"""
function add_cuts!(
    env::DisjunctiveBendersEnv,
    dual_info::Tuple,
    cut_strategy::DisjunctiveCut,
    cut_type::CutType
)
    add_cuts_with_strategy!(
        env,
        dual_info,
        cut_strategy,
        cut_type == K_CUTS ? (k, v) -> k : (k, v) -> v,  # primary cuts
        cut_type == K_CUTS ? (k, v) -> v : (k, v) -> k   # secondary cuts
    )
end


# classical cut
"""
    build_dcglp_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, 
                    coeff_info::Tuple{Float64, Vector{Float64}, Float64}, cut_strategy::AbstractCutStrategy)

Build DCGLP specific cuts (k, v, and master cuts) based on the coefficients and strategy.
"""
function build_dcglp_cuts(
    env::DisjunctiveBendersEnv,
    coeff_info::Tuple{Float64, Vector{Float64}, Float64},
    cut_strategy::ClassicalCut
)
    coefficients_t, coefficients_x, constant_term = coeff_info
    
    # Build DCGLP cuts
    cut_k = @expression(env.dcglp.model, constant_term*env.dcglp.model[:k₀] + sum(coefficients_x .* env.dcglp.model[:kₓ]) + coefficients_t * env.dcglp.model[:kₜ])
    cut_v = @expression(env.dcglp.model, constant_term*env.dcglp.model[:v₀] + sum(coefficients_x .* env.dcglp.model[:vₓ]) + coefficients_t * env.dcglp.model[:vₜ])
    cut_master = @expression(env.master.model, constant_term + sum(coefficients_x .* env.master.model[:x]) + coefficients_t * env.master.model[:t])
    
    return cut_k, cut_v, cut_master
end

