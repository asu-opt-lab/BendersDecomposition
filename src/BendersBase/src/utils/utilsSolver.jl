export assign_attributes!, register_solver!, get_available_solvers
export get_cplex_optimizer, get_cplex_lp_optimizer
export get_gurobi_optimizer, get_gurobi_lp_optimizer

# Registry for available solvers (populated by extensions)
const AVAILABLE_SOLVERS = Dict{Symbol, Any}()

"""
    register_solver!(name::Symbol, optimizer_factory)

Register a solver as available. Called by extensions when they load.
"""
function register_solver!(name::Symbol, optimizer_factory)
    AVAILABLE_SOLVERS[name] = optimizer_factory
    @debug "Registered solver: $name"
end

"""
    get_available_solvers()

Return a list of available solver names.
"""
function get_available_solvers()
    return collect(keys(AVAILABLE_SOLVERS))
end

"""
    is_solver_available(name::Symbol)

Check if a solver is available.
"""
function is_solver_available(name::Symbol)
    return haskey(AVAILABLE_SOLVERS, name)
end

"""
    get_solver(name::Symbol)

Get the optimizer factory for a registered solver.
"""
function get_solver(name::Symbol)
    if !haskey(AVAILABLE_SOLVERS, name)
        available = join(get_available_solvers(), ", ")
        error("Solver :$name is not available. Available solvers: [$available]. " *
              "Make sure to load the solver package (e.g., `using CPLEX` or `using Gurobi`).")
    end
    return AVAILABLE_SOLVERS[name]
end

"""
    assign_attributes!(model::Model, config::Dict{String,Any})

Assign optimizer and attributes to a JuMP model based on configuration.

# Arguments
- `model::Model`: The JuMP model to configure
- `config::Dict{String,Any}`: Configuration dictionary with:
  - `"solver"`: Solver name ("CPLEX", "Gurobi", "HiGHS", etc.)
  - Other keys: Solver-specific parameter names and values

# Example
```julia
config = Dict(
    "solver" => "CPLEX",
    "CPXPARAM_Threads" => 4,
    "CPX_PARAM_EPGAP" => 1e-6
)
assign_attributes!(model, config)
```
"""
function assign_attributes!(model::Model, config::Dict{String,Any})
    solver_name = Symbol(config["solver"])
    
    if !is_solver_available(solver_name)
        available = join(get_available_solvers(), ", ")
        error("Solver $(config["solver"]) is not available. Available solvers: [$available]. " *
              "Make sure to load the solver package (e.g., `using CPLEX` or `using Gurobi`).")
    end
    
    optimizer_factory = get_solver(solver_name)
    set_optimizer(model, optimizer_factory)

    # Set solver attributes from config
    for (param, value) in config
        if param != "solver"
            set_optimizer_attribute(model, param, value)
        end
    end

    set_optimizer_attribute(model, MOI.Silent(), true)
end

# Solver availability flags - set by extensions when loaded
const CPLEX_AVAILABLE = Ref(false)
const GUROBI_AVAILABLE = Ref(false)

# Optimizer factory storage - populated by extensions
const CPLEX_OPTIMIZER_FACTORY = Ref{Any}(nothing)
const CPLEX_LP_OPTIMIZER_FACTORY = Ref{Any}(nothing)
const GUROBI_OPTIMIZER_FACTORY = Ref{Any}(nothing)
const GUROBI_LP_OPTIMIZER_FACTORY = Ref{Any}(nothing)

"""
    get_cplex_optimizer(; kwargs...)

Create a CPLEX optimizer. Requires CPLEX.jl to be loaded.
"""
function get_cplex_optimizer(; kwargs...)
    if !CPLEX_AVAILABLE[]
        error("CPLEX is not available. Please add CPLEX to your project and load it with `using CPLEX`.")
    end
    return CPLEX_OPTIMIZER_FACTORY[](; kwargs...)
end

"""
    get_cplex_lp_optimizer(; kwargs...)

Create a CPLEX optimizer optimized for LP. Requires CPLEX.jl to be loaded.
"""
function get_cplex_lp_optimizer(; kwargs...)
    if !CPLEX_AVAILABLE[]
        error("CPLEX is not available. Please add CPLEX to your project and load it with `using CPLEX`.")
    end
    return CPLEX_LP_OPTIMIZER_FACTORY[](; kwargs...)
end

"""
    get_gurobi_optimizer(; kwargs...)

Create a Gurobi optimizer. Requires Gurobi.jl to be loaded.
"""
function get_gurobi_optimizer(; kwargs...)
    if !GUROBI_AVAILABLE[]
        error("Gurobi is not available. Please add Gurobi to your project and load it with `using Gurobi`.")
    end
    return GUROBI_OPTIMIZER_FACTORY[](; kwargs...)
end

"""
    get_gurobi_lp_optimizer(; kwargs...)

Create a Gurobi optimizer optimized for LP. Requires Gurobi.jl to be loaded.
"""
function get_gurobi_lp_optimizer(; kwargs...)
    if !GUROBI_AVAILABLE[]
        error("Gurobi is not available. Please add Gurobi to your project and load it with `using Gurobi`.")
    end
    return GUROBI_LP_OPTIMIZER_FACTORY[](; kwargs...)
end
