export infeasibility_report

include("utilsSolver.jl")
include("utilsCuts.jl")
include("utilsLoop.jl")
include("utilsBnB.jl")
include("utilsInterface.jl")

"""
    infeasibility_report(master::AbstractMaster, x_opt, t_opt)

Generate and display an infeasibility/consistency check for a candidate solution
to a master model.

This function takes a proposed solution `(x_opt, t_opt)` for the master problem,
loads it into the JuMP model, and evaluates:

  * primal feasibility of all constraints,
  * the objective value at the candidate point,
  * the objective value after fixing all variables and re-solving the model.

The procedure is typically used for debugging to verify
whether the candidate solution is feasible for the master problem.

# Arguments
- `master::AbstractMaster`  
  A container describing the master problem. It must provide:
  - `master.model` — the underlying JuMP model,
  - `master.x` — a collection of `VariableRef`s for the master x-variables,
  - `master.t` — a collection of `VariableRef`s for the t-variables,
  - `master.dim_x` — dimension of `x`,
  - `master.dim_t` — dimension of `t`,
  - `master.c_x`, `master.c_t` — cost vectors for objective evaluation.

- `x_opt::AbstractVector{<:Real}`  
  Candidate values for the x-variables.

- `t_opt::AbstractVector{<:Real}`  
  Candidate values for the t-variables (summed internally to a scalar).

# Behavior
1. Converts the candidate values into a dictionary `opt_sol::Dict{VariableRef,Float64}`.  
2. Prints:
   - a primal feasibility report (`primal_feasibility_report`),  
   - the objective value `cₓᵀ x_opt + cₜᵀ t_opt`.  
3. Fixes all master variables to the candidate values and re-solves the model.  
4. Prints the resulting objective value.

# Returns
Nothing.  
The function is used for logging and diagnostic purposes.

# Side Effects
- Mutates the JuMP model by fixing all variables (`fix(...; force=true)`).
- Solves the model once using `optimize!`.
- Produces logging output via `@info`.

# Example
```julia
# Suppose `master` is an AbstractMaster with dim_x = 3, dim_t = 1
x_opt = [1.0, 0.5, 2.0]
t_opt = [0.3]

infeasibility_report(master, x_opt, t_opt)
"""
function infeasibility_report(master::AbstractMaster, x_opt, t_opt)
    t_opt_ = [sum(t_opt)]
    
    opt_sol = Dict{VariableRef, Float64}()
    for i = 1:master.dim_x
        opt_sol[master.x[i]] = x_opt[i]
    end
    for i = 1:master.dim_t
        opt_sol[master.t[i]] = t_opt_[i]
    end
    
    @info primal_feasibility_report(master.model, opt_sol)
    @info master.c_x' * x_opt + master.c_t' * t_opt_
    
    for v in keys(opt_sol)
        fix(v, opt_sol[v]; force=true)
    end
    optimize!(master.model)
    @info objective_value(master.model)
end