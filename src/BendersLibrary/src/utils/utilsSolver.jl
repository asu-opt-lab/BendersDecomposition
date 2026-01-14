# Solver utilities for BendersLibrary
# Note: Solver-specific implementations are provided by extensions (ext/)
# This file contains solver-agnostic utilities

# The actual optimizer factory functions (get_cplex_optimizer, get_gurobi_optimizer, etc.)
# are defined in the main BendersLibrary.jl file as stubs,
# and are overridden by extensions when CPLEX.jl or Gurobi.jl are loaded.
