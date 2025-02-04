# # Beginner Tutorial

# ## Why should you use BendersDecomposition.jl?

# For the beginners, BendersDecomposition.jl is a good choice for solving optimization problems using Benders Decomposition. It provides an easy-to-use interface for building the model and directly solve the problem using some classical algorithms.
# That means you can solve your problem without worrying about the details of the algorithm.

# We provide some classical algorithms for solving the Benders Decomposition, such as the classical Benders cuts [paper link], the unified Benders cuts [paper link], the Disjunctive Benders cuts [paper link], etc.

# Also, for beginners, we provide both iterative and branch-and-cut algorithms for solving the problems.

# ## @benders_decomposition

# The @benders_decomposition macro is the core of BendersDecomposition.jl. It is used to build the model and solve the problem.

# ## Simple Example

# Here is an example of solving a simple CFLP problem using BendersDecomposition.jl.

# ```julia
# using BendersDecomposition
# using JuMP

# @benders_decomposition standard_env begin
#     @master_problem begin
#         @varible
#         @objective
#         @constraint
#     end
#     @sub_problem begin
#         @variable
#         @objective
#         @constraint
#         @constraint
#     end
# end
# solve!(standard_env)
# ```
