using Test
using JuMP
using Gurobi, CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition

# global_logger(ConsoleLogger(stderr, Logging.Warn))

@testset "CFLP Disjunctive System Tests" begin
    # solver = :Gurobi
    solver = "CPLEX"
    
    # Test on a few representative instances
    for i in [1:66;68:71]
    # for i in 29:66
    # for i in [25]
        @testset "Instance: p$(i)" begin
            # Load CFLP data
            data = read_cflp_benchmark_data("p$(i)")
            
            # Solve using standard MIP model for comparison
            milp = create_milp(data)
            set_optimizer(milp.model, CPLEX.Optimizer)
            optimize!(milp.model)
            mip_objective = objective_value(milp.model)
            
            loop_strategy = Callback()
            disjunctive_system = DisjunctiveCut(ClassicalCut(), L1Norm(), PureDisjunctiveCut(), true, true, true,true)
            # disjunctive_system = DisjunctiveCut(KnapsackCut(), L1Norm(), PureDisjunctiveCut(), true, true, true,true)
            # disjunctive_system = DisjunctiveCut(KnapsackCut(), LInfNorm(), PureDisjunctiveCut(), true, false,false,false)
            
            params = BendersParams(
                60.0,
                1e-5, # *100 already
                solver,
                Dict("solver" => solver),
                Dict("solver" => solver),
                Dict("solver" => solver),
                # Dict(:solver => :Gurobi),
                true
                # false
            )
            benders_UB = Dict()
            benders_LB = Dict()

            obj_value, _ = run_Benders(data, loop_strategy, disjunctive_system, params)
            @test isapprox(mip_objective, obj_value, atol=0.1)

        end
    end
end
