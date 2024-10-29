

struct CallbackBenders <: AbstractBendersAlgorithm
    master_env::AbstractMasterEnv
    sub_env::AbstractSubEnv
    callback::Function
    params::BendersParams
end

function initialize!(algo::CallbackBenders)
    
end

function execute_algorithm(algo::CallbackBenders)
    
end

function finalize!(algo::CallbackBenders, result)
    
end
