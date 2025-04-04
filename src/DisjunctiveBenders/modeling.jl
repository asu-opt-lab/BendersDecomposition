export DisjunctiveBendersEnv
export @disjunctive_benders_decomposition

mutable struct DisjunctiveBendersEnv <: AbstractBendersEnv
    master::AbstractMasterProblem
    sub::AbstractSubProblem
    dcglp::DCGLP
    DisjunctiveBendersEnv() = new()
    DisjunctiveBendersEnv(master::AbstractMasterProblem, sub::AbstractSubProblem, dcglp::DCGLP) = new(master, sub, dcglp)
end


"""
    @disjunctive_benders_decomposition(env_expr, body)

Create a disjunctive Benders decomposition environment.

# Arguments
- `env_expr`: The variable name to store the environment
- `body`: The code block containing master and sub problem definitions and optional configuration

# Example
```julia
# Using default configuration
@disjunctive_benders_decomposition env begin
    @master_problem master begin
        # master problem definition
    end
    @sub_problem sub begin
        # sub problem definition
    end
end

# Using custom configuration
@disjunctive_benders_decomposition env begin
    @master_problem master begin
        # master problem definition
    end
    @sub_problem sub begin
        # sub problem definition
    end
    @config DisjunctiveCut(
        ClassicalCut(),
        L2Norm(),
        StrengthenedDisjunctiveCut(),
        true,
        true,
        true,
        true
    )
end
```
"""
macro disjunctive_benders_decomposition(env_expr, body)
    @assert body.head == :block "Disjunctive Benders decomposition definition must be a code block"
    
    # Create a new variable to store the environment
    env_sym = gensym("env")
    
    # Default configuration
    default_config = quote
        DisjunctiveCut(
            ClassicalCut(),
            L1Norm(),
            PureDisjunctiveCut(),
            false,
            true,
            false,
            false
        )
    end
    
    # Extract configuration if present
    config = default_config
    new_body = Expr(:block)
    
    for expr in body.args
        if expr isa Expr && expr.head == :macrocall && expr.args[1] == Symbol("@config")
            config = expr.args[3]
        elseif expr isa Expr && expr.head == :macrocall
            if expr.args[1] == Symbol("@master_problem") || expr.args[1] == Symbol("@sub_problem")
                push!(new_body.args, Expr(:macrocall, expr.args[1], expr.args[2], env_sym, expr.args[3:end]...))
            else
                push!(new_body.args, expr)
            end
        else
            push!(new_body.args, expr)
        end
    end
    
    return quote
        $(esc(env_expr)) = let
            $(esc(env_sym)) = DisjunctiveBendersEnv()
            $(esc(new_body))
            disjunction_system = $(esc(config))
            $(esc(env_sym)).dcglp = create_dcglp($(esc(env_sym)).master, disjunction_system)
            $(esc(env_sym))
        end
    end
end

