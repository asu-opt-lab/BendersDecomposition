export DisjunctiveOracle

mutable struct DisjunctiveOracle <: AbstractDisjunctiveOracle
    dcglp::Model
    
    # setting
    typical_oracles::Vector{AbstractTypicalOracle}
    norm::AbstractNorm
    split_index_selection_rule::SplitIndexSelectionRule
    disjunctive_cut_append_rule::DisjunctiveCutsAppendRule
    strengthened::Bool
    add_benders_cuts_to_master::Bool
    fraction_of_benders_cuts_to_master::Float64
    reuse_dcglp::Bool
    verbose::Bool

    # log
    disjunctiveCutsByIndex::Vector{Vector{Hyperplane}}
    disjunctiveCuts::Vector{Hyperplane}
    splits::Vector{Tuple{SparseVector{Float64, Int}, Float64}}

    function DisjunctiveOracle(data,  
                            typical_oracles::Vector{T}, 
                            norm::AbstractNorm, 
                            split_index_selection_rule::SplitIndexSelectionRule, 
                            disjunctive_cut_append_rule::DisjunctiveCutsAppendRule; 
                            strengthened::Bool=true, 
                            add_benders_cuts_to_master::Bool=true, 
                            fraction_of_benders_cuts_to_master::Float64 = 1.0, 
                            reuse_dcglp::Bool=true, 
                            verbose::Bool=true) where {T<:AbstractTypicalOracle}
        @debug "Building disjunctive oracle"
        dcglp = Model()
        # Define variables
        @variable(dcglp, tau)
        @variable(dcglp, omega_0[1:2]) # 1 for kappa; 2 for nu
        @variable(dcglp, omega_x[1:2,1:data.dim_x])
        @variable(dcglp, omega_t[1:2,1:data.dim_t])
        @variable(dcglp, sx[1:data.dim_x])
        @variable(dcglp, st[1:data.dim_t])
        
        # Set objective
        @objective(dcglp, Min, tau)

        # Add constraints
        @constraint(dcglp, [i=1:2], omega_t[i,:] .>= -1e6 * omega_0[i])
        @constraint(dcglp, coneta[i in 1:2, j in 1:data.dim_x], 0 >= -omega_0[i] + omega_x[i,j]) 
        @constraint(dcglp, condelta[i in 1:2, j in 1:data.dim_x], 0 >= -omega_x[i,j])
        @constraint(dcglp, conineq[i in 1:2], omega_0[i] >= 0)

        # Add gamma constraints
        @constraint(dcglp, con0, omega_0[1] + omega_0[2] == 1)
        @constraint(dcglp, conx, omega_x[1,:] + omega_x[2,:] - sx .== 0)
        @constraint(dcglp, cont[j=1:data.dim_t], omega_t[1,j] + omega_t[2,j] - st[j] == 0) # must be in this form to recognize it as a vector

        add_normalization_constraint(data, dcglp, norm)

        if disjunctive_cut_append_rule == DisjunctiveCutsSmallerIndices()
            if !(typeof(split_index_selection_rule) <: SimpleSplit)
                throw(AlgorithmException("$(typeof(disjunctive_cut_append_rule)) can only be paired with SimpleSplit"))
            end
        end
        
        disjunctiveCutsByIndex = [Vector{Hyperplane}() for i=1:data.dim_x]
        splits = Vector{Tuple{SparseVector{Float64, Int}, Float64}}()

        new(dcglp, typical_oracles, norm, split_index_selection_rule, disjunctive_cut_append_rule, strengthened, add_benders_cuts_to_master, fraction_of_benders_cuts_to_master, reuse_dcglp, verbose, disjunctiveCutsByIndex, Vector{Hyperplane}(), splits)
    end
end

function generate_cuts(oracle::DisjunctiveOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol = 1e-6, time_limit = 3600)

    tic = time()
    
    push!(oracle.splits, select_disjunctive_inequality(x_value, oracle.split_index_selection_rule))
    index = get_split_index(oracle)
    
    if get_sec_remaining(tic, time_limit) <= 0.0
        throw(TimeLimitException("Time limit reached during cut generation"))
    end

    replace_disjunctive_inequality!(oracle)
    
    # delete benders cuts previously added when not reusing dcglp
    if !oracle.reuse_dcglp
        if haskey(oracle.dcglp, :con_benders)
            delete(oracle.dcglp, oracle.dcglp[:con_benders]) 
            unregister(oracle.dcglp, :con_benders)
        end
    end

    # add previously found disjunctive cuts based on a user-given append rule
    add_disjunctive_cuts!(oracle, oracle.disjunctive_cut_append_rule)

    if get_sec_remaining(tic, time_limit) <= 0.0
        throw(TimeLimitException("Time limit reached during cut generation"))
    end

    set_normalized_rhs.(oracle.dcglp[:conx], x_value)
    set_normalized_rhs.(oracle.dcglp[:cont], t_value)

    return solve_dcglp!(oracle, x_value, t_value; time_limit = time_limit)
end

"""
prototypes for user-customizable functions for DisjunctiveOracle
"""
function add_normalization_constraint(data::Data, dcglp::Model, norm::AbstractNorm)
    throw(UndefError("update add_normalization_constraint for $(typeof(norm))"))
    # should add a normalization constraint to dcglp. 
end

function select_disjunctive_inequality(x_value::Vector{Float64}, split_selection_rule::SplitIndexSelectionRule; zero_tol = 1e-2)    
    throw(UndefError("update select_disjunctive_inequality for $(typeof(split_selection_rule))"))
    # should return a split: phi, phi_0
end

function add_disjunctive_cuts!(oracle::DisjunctiveOracle, rule::DisjunctiveCutsAppendRule)
    throw(UndefError("update add_disjunctive_cuts! for $(typeof(rule))"))
end

include("oracleDisjunctiveInterface.jl")

"""
utility functions for DisjunctiveOracle
"""
function get_split_index(oracle::DisjunctiveOracle)
    if !(typeof(oracle.split_index_selection_rule) <: SimpleSplit)
        throw(AlgorithmException("get_split_index is only valid for SimpleSplit"))
    end
    return findfirst(x -> x > 0.5, oracle.splits[end][1])
end

function replace_disjunctive_inequality!(oracle::DisjunctiveOracle)
    dcglp = oracle.dcglp
    phi = oracle.splits[end][1]
    phi_0 = oracle.splits[end][2]
    
    if haskey(dcglp, :con_split_kappa)
        delete(dcglp, dcglp[:con_split_kappa]) 
        unregister(dcglp, :con_split_kappa)
    end
    if haskey(dcglp, :con_split_nu)
        delete(dcglp, dcglp[:con_split_nu]) 
        unregister(dcglp, :con_split_nu)
    end
        
    # Add new constraints
    @constraint(dcglp, con_split_kappa, 0 >= dcglp[:omega_0][1]*(phi_0+1) - phi' * dcglp[:omega_x][1,:])
    @constraint(dcglp, con_split_nu, 0 >= -dcglp[:omega_0][2]*phi_0 + phi' * dcglp[:omega_x][2,:])
end

function solve_dcglp!(oracle::DisjunctiveOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; zero_tol = 1e-9, time_limit = time_limit)

    # setting: halt_limit = 3, 
    log = DcglpLog()
    state = DcglpState()
    dcglp = oracle.dcglp

    # 1 for k; 2 for nu
    typical_oracles = oracle.typical_oracles

    # cuts for master
    hyperplanes = Vector{Hyperplane}()

    prev_lb = -Inf

    while true
        state.iteration += 1

        state.master_time = @elapsed begin
            set_time_limit_sec(dcglp, get_sec_remaining(log.start_time, time_limit))
            optimize!(dcglp)
            if is_solved_and_feasible(dcglp; allow_local = false, dual = true)
                for i=1:2
                    omega_value[:x][i] = value.(dcglp[:omega_x][i,:])
                    omega_value[:t][i] = value.(dcglp[:omega_t][i,:])
                    omega_value[:z][i] = value(dcglp[:omega_0][i])
                end
                other_values = (tau = value(dcglp[:tau]), sx = value.(dcglp[:sx]))
                state.LB = other_values.tau
            elseif termination_status(dcglp) == ALMOST_INFEASIBLE
                @warn "dcglp master termination status: $(termination_status(dcglp)); the problem is infeasible or dcglp encountered numerical issue, yielding the typical Benders cut"
                return generate_cuts(typical_oracles[1], x_value, t_value; time_limit = get_sec_remaining(log.start_time, time_limit))
            elseif termination_status(env.master.model) == TIME_LIMIT
                throw(TimeLimitException("Time limit reached during dcglp solving"))
            else
                throw(UnexpectedModelStatusException("dcglp master: $(termination_status(dcglp))"))
                # if infeasible, then the problem is infeasible
            end
        end
        # log.master_time += state.master_time
        
        benders_cuts = Dict(1 => Vector{AffExpr}(), 2 => Vector{AffExpr}())
        # my_lock = Threads.ReentrantLock()
        # Threads.@threads for i in 1:2
        for i in 1:2
            # Threads.lock(my_lock) do
                state.oracle_times[i] = @elapsed begin
                    if omega_value[:z][i] >= zero_tol
                        state.is_in_L[i], hyperplanes_a, state.f_x[i] = generate_cuts(typical_oracles[i], omega_value[:x][i] / omega_value[:z][i], omega_value[:t][i] / omega_value[:z][i], time_limit = get_sec_remaining(log.start_time, time_limit))
                        # adjust the tolerance with respect to dcglp: (sum(state.sub_obj_vals[i]) - sum(t_value)) * omega_value[:z][i] < zero_tol
                        if !state.is_in_L[i]
                            for k = 1:2 # add to both kappa and nu systems
                                append!(benders_cuts[i], hyperplanes_to_expression(dcglp, hyperplanes_a, dcglp[:omega_x][k,:], dcglp[:omega_t][k,:], dcglp[:omega_0][k]))
                            end
                            if oracle.add_benders_cuts_to_master
                                append!(hyperplanes, select_top_fraction(hyperplanes_a, h -> evaluate_violation(h, x_value, t_value), oracle.fraction_of_benders_cuts_to_master))
                            end
                        end
                    else
                        state.is_in_L[i] = true
                        state.sub_obj_vals[i] = zeros(length(t_value))
                    end
                end
                # log.oracle_times[i] += state.oracle_times[i]
            # end
        end
    
        if state.sub_obj_vals[1] !== nothing && state.sub_obj_vals[2] !== nothing
            update_upper_bound_and_gap!(state, omega_value, other_values, t_value, oracle.norm)
        end
        
        push!(log.iterations, state)

        oracle.verbose && print_dcglp_iteration_info(state)

        check_lb_improvement!(state, prev_lb; zero_tol = zero_tol)

        is_terminated(state, log, time_limit) && break
        
        add_constraints(dcglp, :con_benders, [benders_cuts[1]; benders_cuts[2]]) 

        prev_lb = state.LB
    end

    if state.LB >= zero_tol
        gamma_x, gamma_t, gamma_0 = oracle.strengthened ? generate_strengthened_disjunctive_cuts(oracle.dcglp) : generate_disjunctive_cut(oracle.dcglp)

        h = Hyperplane(gamma_x, gamma_t, gamma_0)
        push!(hyperplanes, h)
        
        if typeof(oracle.split_index_selection_rule) <: SimpleSplit
            index = get_split_index(oracle)
            push!(oracle.disjunctiveCutsByIndex[index], h)
        end
        push!(oracle.disjunctiveCuts, h)
        
        if oracle.disjunctive_cut_append_rule == AllDisjunctiveCuts()
            d_cuts = Vector{AffExpr}()
            for k = 1:2 # add to both kappa and nu systems
                append!(d_cuts, hyperplanes_to_expression(dcglp, [h], dcglp[:omega_x][k,:], dcglp[:omega_t][k,:], dcglp[:omega_0][k]))
            end
            add_constraints(dcglp, :con_disjunctive, d_cuts) 
        end
        
        return false, hyperplanes, fill(Inf, length(t_value))
    else
        return generate_cuts(typical_oracles[1], x_value, t_value; time_limit = get_sec_remaining(log.start_time, time_limit))
    end
    # statistics_of_disjunctive_cuts(env)
end

function generate_disjunctive_cut(dcglp::Model)
    gamma_x = dual.(dcglp[:conx])
    gamma_t = dual.(dcglp[:cont])
    gamma_0 = dual(dcglp[:con0])
    
    return gamma_x, gamma_t, gamma_0
end

function generate_strengthened_disjunctive_cuts(dcglp::Model; zero_tol = 1e-5)
    
    σ₁ = dual(dcglp[:con_split_kappa])
    σ₂ = dual(dcglp[:con_split_nu])
    gamma_x = dual.(dcglp[:conx])
    gamma_t = dual.(dcglp[:cont])
    gamma_0 = dual(dcglp[:con0])

    # println("DCGLP Sigma Values: [σ₁: $σ₁, σ₂: $σ₂]")

    a₁ = -gamma_x .- dual.(dcglp[:condelta][1])
    a₂ = -gamma_x .- dual.(dcglp[:condelta][2])
    σ_sum = σ₂ + σ₁
    if σ_sum >= zero_tol
        m = (a₁ .- a₂) / σ_sum
        m_lb = floor.(m)
        m_ub = ceil.(m)
        gamma_x = -min.(a₁-σ₁*m_lb, a₂+σ₂*m_ub)
    end
    
    return gamma_x, gamma_t, gamma_0
end

"""
utility functions for Dcglp
"""
mutable struct DcglpState <: AbstractDcglpState
    master_time::Float64
    oracle_times::Vector{Float64}
    omega_values::Dict{Symbol,Vector}
    other_values::Dict{Symbol,Any}
    f_x::Vector{Vector{Float64}}
    omega_t_::Vector{Vector{Float64}}
    is_in_L::Vector{Bool}
    LB::Float64
    UB::Float64
    gap::Float64

    # Constructor with default values
    function DcglpState() 
        new(0.0, 
            [0.0; 0.0], 
            Dict(:x => Vector{Vector{Float64}}(undef, 2), 
                 :t => Vector{Vector{Float64}}(undef, 2), 
                 :z => Vector{Float64}(undef, 2)), 
            Dict(:tau => -Inf, :sx => Vector{Float64}()), 
            Vector{Vector{Float64}}(undef, 2), 
            Vector{Vector{Float64}}(undef, 2), 
            [false; false], 
            -Inf, 
            Inf, 
            100.0)
    end
end

mutable struct DcglpLog <: AbstractDcglpLog
    n_iter::Int
    iterations::Vector{DcglpState}
    start_time::Float64
    consecutive_no_improvement::Int
    # master_time::Float64 # do we need this?
    # oracle_times::Vector{Float64}
    
    function DcglpLog()
        new(0, Vector{DcglpState}(), time(), 0)
    end
end

# # the followings are for loop: 
# function update_upper_bound_and_gap!(state::AbstractLoopState, omega_values, other_values, t_value, f::function)
# end

function update_upper_bound_and_gap!(state::DcglpState, omega_values, other_values, t_value, norm::AbstractNorm)
    throw(UndefError("update update_upper_bound_and_gap! for $(typeof(norm))"))
end

function update_upper_bound_and_gap!(state::DcglpState, omega_values, other_values, t_value, norm::LpNorm)
    state.omega_values
    state.other_values
    for i=1:2 
        state.omega_t_[i] = state.is_in_L[i] ? state.omega_values[:t][i] : state.f_x[i] * state.omega_values[:z][i]
    end
    st_ = state.omega_t_[1] .+ state.omega_t_[2] .- t_value
    evaluation = f(st_, state.t_value)
    (t1, t2) -> LinearAlgebra.norm([other_values.sx; t1 .+ t2 .- t_value], norm.p)
    state.UB = min(state.UB, LinearAlgebra.norm([other_values.sx; st_], norm.p))
    state.gap = (state.UB - state.LB) / abs(state.UB) * 100
end

"""
Print iteration information if verbose mode is on
"""
function print_iteration_info(state::DcglpState, log::DcglpLog)
    @printf("   Iter: %4d | LB: %8.4f | UB: %8.4f | Gap: %6.2f%% | UB_k: %8.2f | UB_v: %8.2f | Master time: %6.2f | Sub_k time: %6.2f | Sub_v time: %6.2f \n",
           log.n_iter, state.LB, state.UB, state.gap, sum(state.omega_t_[1]), sum(state.omega_t_[2]), state.master_time, state.oracle_times[1], state.oracle_times[2])
end

"""
Check termination criteria
"""
function is_terminated(state::DcglpState, log::DcglpLog, params::BendersParams; time_limit::Float64)
    return state.is_in_L[1] && state.is_in_L[2] || log.consecutive_no_improvement >= param.dcglp.halt_limit || state.gap <= params.dcglp.gap_tolerance  || get_sec_remaining(log.start_time, time_limit) <= 0.0 || time() - log.start_time >= param.dcglp.time_limit || log.n_iter >= param.dcglp.iter_limit
end








