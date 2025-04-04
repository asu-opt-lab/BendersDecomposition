using Test
using JuMP
using CPLEX
using Printf
using DataFrames
using Logging
using BendersOracle

#################
# TODO:
#################
@testset "CFLP Sequential Benders Tests" begin
    solver = "CPLEX"
    # solver = :Gurobi
    # instances = [1:66; 68:71]
    instances = [20]
    for i in instances
        @testset "Instance: p$i" begin
            # Load data
            data = read_cflp_benchmark_data("p$i")
            # data = read_GK_data("f100-c100-r5-1")
            # Create and solve MIP reference model
            milp = create_milp(data)
            set_optimizer(milp.model, CPLEX.Optimizer)
            optimize!(milp.model)
            mip_obj = objective_value(milp.model)
            @info mip_obj

            cut_strategy = ClassicalCut()
            loop_strategy = GenericCallback(100, 1e-6, lazy_callback, nothing, true)
            master = create_master_problem(data, cut_strategy)
            sub = create_sub_problem(data, cut_strategy)
            env = BendersEnv(master, sub)
            solve!(env, loop_strategy, cut_strategy)
        end
    end
end