import Dualization
import SCS
import Ipopt

export DCGLP

mutable struct UFLPDCGLPEnv <: AbstractDCGLPEnv
    model::Model
    ifsolved::Bool
    masterconπpoints1
    masterconπpoints2
    conπpoints1
    conπpoints2
end

function DCGLP(sub_env::AbstractSubEnv, a::Vector{Int}, b::Int, norm_type::AbstractNormType; solver::Symbol=:Gurobi)
    @error "wrong type set"
end

# function DCGLP(sub_env::UFLPSplitSubEnv, a::Vector{Int}, b::Int, norm_type::StandardNorm; solver::Symbol=:Gurobi) 

#     data = sub_env.data
#     if solver == :CPLEX
#         model = Model(CPLEX.Optimizer)
#         # set_optimizer_attribute(model, "CPX_PARAM_LPMETHOD", CPX_ALG_DUAL)
#     else
#         model = Model(Gurobi.Optimizer)
#         # set_optimizer_attribute(model, "InfUnbdInfo", 1)
#     end
#     model =  Model(CPLEX.Optimizer)

#     set_optimizer_attribute(model, MOI.Silent(),true)

#     # pre
#     N = data.n_facilities
#     M = data.n_customers
    
#     # Variables
#     @variable(model, τ)
#     @variable(model, k₀>=0)
#     @variable(model, kₓ[1:N])
#     @variable(model, kₜ)
#     @variable(model, v₀>=0)
#     @variable(model, vₓ[1:N])
#     @variable(model, vₜ)

#     # Objective
#     @objective(model, Min, τ)


#     # Constraints
#     total_demands = sum(data.demands)

#     @constraint(model, consigma1, τ >= k₀*(b+1) - a'kₓ) 
#     @constraint(model, coneta1[j in 1:N], τ >= -k₀ + kₓ[j]) 
#     @constraint(model, consigma2, τ >= -v₀*b + a'vₓ ) 
#     @constraint(model, coneta2[j in 1:N], τ >= -v₀ + vₓ[j])

#     @constraint(model, conv1[j in 1:N], τ >= -kₓ[j])
#     @constraint(model, conw1, τ >= -sum(data.capacities[j]*kₓ[j] for j in 1:N) + total_demands*k₀)

#     @constraint(model, conv2[j in 1:N], τ >= -vₓ[j])
#     @constraint(model, conw2, τ >= -sum(data.capacities[j]*vₓ[j] for j in 1:N) + total_demands*v₀)

#     @constraint(model, con0, k₀ + v₀ == 1)
#     @constraint(model, conx[i=1:N], kₓ[i] .+ vₓ[i] == 0)  #-x̂
#     @constraint(model, cont, kₜ + vₜ == 0) #-t̂
    
#     for i in eachindex(sub_env.split_info.γ₀s)
#         @constraint(model, sub_env.split_info.γ₀s[i]*k₀ + sub_env.split_info.γₓs[i]'kₓ + sub_env.split_info.γₜs[i]*kₜ <= 0)
#         @constraint(model, sub_env.split_info.γ₀s[i]*v₀ + sub_env.split_info.γₓs[i]'vₓ + sub_env.split_info.γₜs[i]*vₜ <= 0)
#     end

#     return CFLPDCGLPEnv(model, false, [], [], [], [])
# end
    


# function DCGLP(sub_env::CFLPSplitSubEnv, a::Vector{Int}, b::Int, norm_type::GammaNorm; solver::Symbol=:Gurobi)
function DCGLP(sub_env::UFLPSplitSubEnv, a::Vector{Int}, b::Int, norm_type::GammaNorm; solver::Symbol=:CPLEX)
    data = sub_env.data
    if solver == :CPLEX
        model = Model(CPLEX.Optimizer)
        # set_optimizer_attribute(model, "CPX_PARAM_LPMETHOD", CPX_ALG_BARRIER)
    elseif solver == :Gurobi
        model = Model(Gurobi.Optimizer)
        # set_optimizer_attribute(model, "InfUnbdInfo", 1)
    end

    set_optimizer_attribute(model, MOI.Silent(),true)

    # pre
    N = data.n_facilities
    M = data.n_customers
    
    # Variables
    @variable(model, τ)
    @variable(model, k₀>=0)
    @variable(model, kₓ[1:N])
    @variable(model, kₜ)
    @variable(model, v₀>=0)
    @variable(model, vₓ[1:N])
    @variable(model, vₜ)
    @variable(model, sx[i = 1:N])
    @variable(model, st)

    # Objective
    @objective(model, Min, τ)
    
    total_demands = sum(data.demands)

    # Constraints
    @constraint(model, consigma1, 0 >= k₀*(b+1) - a'kₓ) 
    @constraint(model, coneta1[j in 1:N], 0 >= -k₀ + kₓ[j]) 
    @constraint(model, consigma2, 0 >= -v₀*b + a'vₓ) 
    @constraint(model, coneta2[j in 1:N], 0 >= -v₀ + vₓ[j])

    @constraint(model, conv1[j in 1:N], 0 >= -kₓ[j])
    # @constraint(model, conw1, 0 >= -sum(data.capacities[j]*kₓ[j] for j in 1:N) + total_demands*k₀)

    @constraint(model, conv2[j in 1:N], 0 >= -vₓ[j])
    # @constraint(model, conw2, 0 >= -sum(data.capacities[j]*vₓ[j] for j in 1:N) + total_demands*v₀)

    @constraint(model, con0, k₀ + v₀ == 1)
    if norm_type == L1GAMMANORM
        @constraint(model, concone, [τ; sx ; st] in MOI.NormInfinityCone(2 + N))
    elseif norm_type == L2GAMMANORM
        @constraint(model, concone, [τ; sx ; st] in MOI.SecondOrderCone(2 + N))
    elseif norm_type == LINFGAMMANORM
        @constraint(model, concone, [τ; sx ; st] in MOI.NormOneCone(2 + N))
    end
    
    @constraint(model, conx[i=1:N], kₓ[i] .+ vₓ[i] .- sx[i] == 0)  #x̂
    @constraint(model, cont, kₜ + vₜ - st == 0) #t̂
    
    
    for i in eachindex(sub_env.split_info.γ₀s)
        @constraint(model, sub_env.split_info.γ₀s[i]*k₀ + sub_env.split_info.γₓs[i]'kₓ + sub_env.split_info.γₜs[i]*kₜ <= 0)
        @constraint(model, sub_env.split_info.γ₀s[i]*v₀ + sub_env.split_info.γₓs[i]'vₓ + sub_env.split_info.γₜs[i]*vₜ <= 0)
    end

    return UFLPDCGLPEnv(model, false, [], [], [], [])
end

