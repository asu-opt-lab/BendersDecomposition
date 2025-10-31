export SpecializedBendersSeq

mutable struct SpecializedBendersSeq <: AbstractBendersSeq
    data::Data
    master::AbstractMaster
    oracle::DisjunctiveOracle

    param::SpecializedBendersSeqParam

    # result
    obj_value::Float64
    termination_status::TerminationStatus

    function SpecializedBendersSeq(data, master::AbstractMaster, oracle::DisjunctiveOracle; param::SpecializedBendersSeqParam = SpecializedBendersSeqParam()) 
        relax_integrality(master.model)

        oracle.oracle_param.split_index_selection_rule != LargestFractional() && throw(AlgorithmException("SpeicalizedBendersSeq does not admit $(oracle.oracle_param.split_index_selection_rule). Use LargestFractional() instead."))
        oracle.oracle_param.disjunctive_cut_append_rule != DisjunctiveCutsSmallerIndices() && throw(AlgorithmException("SpeicalizedBendersSeq does not admit $(oracle.oracle_param.disjunctive_cut_append_rule). Use DisjunctiveCutsSmallerIndices() instead."))

        new(data, master, oracle, param, Inf, NotSolved())
    end
end

"""
Run SpecializedBendersSeq
"""
function solve!(env::SpecializedBendersSeq) 
    log = BendersSeqLog()
    L_param = BendersSeqParam(; time_limit = env.param.time_limit, gap_tolerance = env.param.lp_gap_tolerance, verbose = env.param.verbose)
    L_env = BendersSeq(env.data, env.master, env.oracle.typical_oracles[1]; param = L_param)

    try
        while true
            state = BendersSeqState()
            state.total_time = @elapsed begin
                ## add all found disjunctive cuts to master
                all_disj_cuts = hyperplanes_to_expression(env.master.model, env.oracle.disjunctiveCuts, env.master.x, env.master.t)    
                @constraint(env.master.model, con_disjunctive, 0.0 .>= all_disj_cuts)

                # Solve linear relaxation
                state.master_time = @elapsed begin
                    solve!(L_env; iter_prefix = " LP")
                    state.LB, state.values[:x], state.values[:t] = JuMP.objective_value(env.master.model), value.(env.master.x), value.(env.master.t)
                end
                @debug state.values[:x]

                # Check termination criteria
                is_terminated(state, log, env.param) && (record_iteration!(log, state); break)

                # generate optimal solution that is vertex of P^{k,j}; essential for numerical stability
                generate_optimal_vertex!(env, L_env, state)
                
                state.oracle_time = @elapsed begin
                    state.is_in_L, hyperplanes, state.f_x = generate_cuts(env.oracle, state.values[:x], state.values[:t]; time_limit = get_sec_remaining(log, env.param), throw_typical_cuts_for_errors = false, include_disjunctive_cuts_to_hyperplanes = false)
                    cuts = !state.is_in_L ? hyperplanes_to_expression(env.master.model, hyperplanes, env.master.x, env.master.t) : []
                end

                state.is_in_L && throw(UnexpectedModelStatusException("SpecializedBendersSeq: τ=0 at fractional point, possibly numerical issue"))
                
                record_iteration!(log, state)

                env.param.verbose && print_iteration_info(state, log)

                # Add generated cuts to master
                @constraint(env.master.model, 0 .>= cuts)
            end
        end
        env.termination_status = Optimal()
        env.obj_value = log.iterations[end].LB
        
        return to_dataframe(log)
    catch e
        if typeof(e) <: TimeLimitException
            env.termination_status = TimeLimit()
            env.obj_value = log.iterations[end].LB
        elseif typeof(e) <: UnexpectedModelStatusException
            env.termination_status = InfeasibleOrNumericalIssue()
            @warn e.msg
        else
            rethrow()  
        end
        if env.param.verbose
            println("Terminated with $(env.termination_status)")
        end
        return to_dataframe(log)
    end
end

function generate_optimal_vertex!(env::SpecializedBendersSeq, L_env::AbstractBendersSeq, state::BendersSeqState)
    dim_x = env.data.dim_x
    
    ## find largest fractional idx
    frac_idx = -1
    for idx in reverse(1:dim_x)
        val = state.values[:x][idx]
        if !isapprox(abs(val - 0.5), 0.5, atol = env.param.integrality_tolerance) 
            frac_idx = idx
            break
        end
    end
    @debug "frac_idx: $frac_idx"

    ## modify master problem
    env.master.model[:fix_x] = @constraint(env.master.model, [i in frac_idx+1:dim_x], env.master.x[i] == round(state.values[:x][i]))
    ## remove and add all disjunctive cuts up to idx
    if haskey(env.master.model, :con_disjunctive)
        delete.(env.master.model, env.master.model[:con_disjunctive]) 
        unregister(env.master.model, :con_disjunctive)
    end
    disj_cuts_idx = reduce(vcat, [hyperplanes_to_expression(env.master.model, env.oracle.disjunctiveCutsByIndex[i], env.master.x, env.master.t) for i = 1:frac_idx])
    @constraint(env.master.model, con_disjunctive, 0.0 .>= disj_cuts_idx)

    ## solve master again
    solve!(L_env; iter_prefix = " LP aux")

    ## compare the solutions; the optimal obj. val. should match
    LB, x_val, t_val = JuMP.objective_value(env.master.model), value.(env.master.x), value.(env.master.t)
    @debug "state.LB = $(state.LB) vs LB = $LB"
    @debug "state.values[:x] = $(state.values[:x]) vs x_val = $x_val"
    !isapprox(state.LB, LB, rtol = 1e-4) && throw(UnexpectedModelStatusException("SpecializedBendersSeq: fail to generate vertex for P^{k,j}, possibly numerical issue"))
    state.LB = LB
    state.values[:x] = x_val
    state.values[:t] = t_val

    ## remove fixing constraints and disjunctive cuts
    if haskey(env.master.model, :fix_x)
        delete.(env.master.model, env.master.model[:fix_x]) 
        unregister(env.master.model, :fix_x)
    end
    if haskey(env.master.model, :con_disjunctive)
        delete.(env.master.model, env.master.model[:con_disjunctive]) 
        unregister(env.master.model, :con_disjunctive)
    end
end