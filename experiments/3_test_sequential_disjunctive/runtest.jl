# To-Do:
# ***** obtain optimal solution, check whether any disjunctive cut cuts off that point; if so, how much.
# add X to dcglp
# streamline the definition for attributes: (1) model vs dcglp; (2) termination parameter for dcglp; (3) verbose for master and oracle directly
# test multiple scenarios -- should work
# add lifting
# continuous master; check split index rule for continuous master
# define BendersSeq and BendersSeqInOut and BendersBnB
# define BendersSeq and BendersSeqInOut and BendersBnB for disjunctive? continuous master?
# plug in other solver, e.g., conic solvers like Mosek 
# add interface for parameters
# param for setting details about log (returned dataframe from solve!: all iterations or last)

# Done:
# when no fractional value exists, randomly choose any index --> done; need to check select_disjunctive_inequality for LP master
# even if disjunctiveOracle does not generate a cut, we should not terminate it --> done; when latest tau is close to zero, return typical Benders cut
# test fat knapsack --> works
# test strengthening technique --> works
# add_benders_cuts_to_master only for a farction of violated ones by (x_value, t_value) --> done
# add disjunctive cut handler for dcglp

# Issues:
# 1. DisjunctiveOracle: Setting CPX_PARAM_EPRHS tightly (e.g., < 1e-6) results in dcglp terminating with ALMOST_INFEASIBLE --> set it tightly (1e-9) and outputs a typical Benders cut when dcglp is ALMOST_INFEASIBLE
# 2. DisjunctiveOracle: Setting zero_tol in solve_dcglp! large (e.g., 1e-4) results in disjunctive Benders with reuse_dcglp = true terminating with incorrect solution --> set it tightly as 1e-9
# 3. DisjunctiveOracle: solve_dcglp! becomes stall since the true violation should be multipled with omega_value[:z], which can be fairly small --> terminate dcglp when LB does not improve for a certain number of iterations
# 4. DisjunctiveOracle: the fat-knapsack-based disjunctive cut may have a sparse gamma_t, so adding only disjunctive cut does not improve lower bound, add_benders_cuts_to_master should be set at true
using Test
using JuMP
using CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition
import BendersDecomposition: generate_cuts


@testset "Sequential Disjunctive Tests" begin
    @info "Running Sequential Disjunctive Tests"
    include("ufl.jl")
    include("cfl.jl")
    include("scfl.jl")
    @info "Sequential Disjunctive Tests completed"
end