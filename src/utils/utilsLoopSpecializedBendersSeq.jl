export SpecializedBendersSeqParam 

mutable struct SpecializedBendersSeqParam <: AbstractBendersSeqParam

    time_limit::Float64
    integrality_tolerance::Float64
    halt_limit::Int
    iter_limit::Int
    verbose::Bool

    function SpecializedBendersSeqParam(; 
                        time_limit::Float64 = 7200.0, 
                        integrality_tolerance::Float64 = 1e-9,
                        halt_limit::Int = Int(1e9), 
                        iter_limit::Int = Int(1e9), 
                        verbose::Bool = true
                        ) 
        
        new(time_limit, integrality_tolerance, halt_limit, iter_limit, verbose)
    end
end

function is_terminated(state::BendersSeqState, log::BendersSeqLog, param::SpecializedBendersSeqParam)
    return all(x -> isapprox(0.5, abs(x-0.5), atol = param.integrality_tolerance), state.values[:x]) || get_sec_remaining(log, param) <= 0.0
end
