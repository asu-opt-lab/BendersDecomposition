module CPLEXWrapper

using CPLEX, MathOptInterface, JuMP
const MOI = MathOptInterface


function submit_local(model::GenericModel, cb::MOI.UserCut, con::ScalarConstraint)
    return submit_local(backend(model), cb, moi_function(con.func), con.set)
end

function submit_local(
    caching_opt::MOI.Utilities.CachingOptimizer,
    sub::MOI.AbstractSubmittable,
    args...,
)
    return submit_local(
        caching_opt.optimizer,
        sub,
        MOI.Utilities.map_indices.(Ref(caching_opt.model_to_optimizer_map), args)...,
    )
end

function submit_local(
    b::MOI.Bridges.AbstractBridgeOptimizer,
    sub::MOI.AbstractSubmittable,
    args...,
)
    return submit_local(b.model, sub, MOI.Bridges.bridged_function.(b, args)...)
end

function submit_local(
    model::CPLEX.Optimizer,
    cb::MOI.UserCut{CPLEX.CallbackContext},
    f::MOI.ScalarAffineFunction{Float64},
    s::Union{
        MOI.LessThan{Float64},
        MOI.GreaterThan{Float64},
        MOI.EqualTo{Float64},
    },
)
    @info "submit"
    if model.callback_state == CPLEX._CB_LAZY
        throw(MOI.InvalidCallbackUsage(MOI.LazyConstraintCallback(), cb))
    elseif model.callback_state == CPLEX._CB_HEURISTIC
        throw(MOI.InvalidCallbackUsage(MOI.HeuristicCallback(), cb))
    elseif !iszero(f.constant)
        throw(
            MOI.ScalarFunctionConstantNotZero{Float64,typeof(f),typeof(s)}(
                f.constant,
            ),
        )
    end
    
    rmatind, rmatval = CPLEX._indices_and_coefficients(model, f)
    sense, rhs = CPLEX._sense_and_rhs(s)
    
    ret = CPLEX.CPXcallbackaddusercuts(
        cb.callback_data,
        Cint(1),
        Cint(length(rmatval)),
        Ref(rhs),
        Ref(sense),
        Ref{Cint}(0),
        rmatind,
        rmatval,
        Ref{Cint}(CPLEX.CPX_USECUT_PURGE),
        Ref{Cint}(1),  
    )
    
    CPLEX._check_ret(model, ret)
    return
end

export submit_local


end # module