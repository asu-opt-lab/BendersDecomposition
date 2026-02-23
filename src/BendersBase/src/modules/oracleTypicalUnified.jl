export UnifiedOracle

mutable struct UnifiedOracle <: AbstractTypicalOracle
    oracle_param::BasicOracleParam
    model::Model

    function UnifiedOracle(oracle_param::BasicOracleParam, model::Model)
        new(oracle_param, model)
    end

    function UnifiedOracle(data::Data;
                           scen_idx::Int=1,
                           solver_param::Dict{String,Any}=Dict(
                               "solver" => "CPLEX",
                               "CPX_PARAM_EPRHS" => 1e-9,
                               "CPX_PARAM_NUMERICALEMPHASIS" => 1,
                               "CPX_PARAM_EPOPT" => 1e-9
                           ),
                           oracle_param::BasicOracleParam=BasicOracleParam())

        @debug "Building unified oracle"
        model = Model()

        # Define coupling variables and constraints
        @variable(model, x[1:data.dim_x])

        assign_attributes!(model, solver_param)

        new(oracle_param, model, gpu_use)
    end

    UnifiedOracle() = new()
end

function generate_cuts(oracle::UnifiedOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    set_time_limit_sec(oracle.model, time_limit)

    set_normalized_rhs.(oracle.model[:fix_x_lb], x_value)
    set_normalized_rhs.(oracle.model[:fix_x_ub], -x_value)
    set_normalized_rhs.(oracle.model[:epigraph], -t_value)

    optimize!(oracle.model)
    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    elseif termination_status(oracle.model) != OPTIMAL
        throw(UnexpectedModelStatusException("UnifiedOracle: $(termination_status(oracle.model)). This is likely a numerical issue."))
    end
    
    a_x = dual.(oracle.model[:fix_x_lb]) .- dual.(oracle.model[:fix_x_ub])
    a_t = [-dual(oracle.model[:epigraph])]
    a_0 = objective_value(oracle.model) - a_x'*x_value + dual(oracle.model[:epigraph])*t_value[1]

    return isapprox(dual_objective_value(oracle.model), 0, atol=oracle.oracle_param.zero_tol) ? (true, [Hyperplane(a_x, a_t, a_0)], [t_value[1]]) : (false, [Hyperplane(a_x, a_t, a_0)], [NaN])
end



