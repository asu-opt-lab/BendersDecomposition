using BendersDecomposition
using Test
using JuMP

benders_param = BendersBnBParam(
    time_limit = 200.0,
    gap_tolerance = 1e-6,
    verbose = true
)
# Solver parameters
mip_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPXPARAM_Threads" => 4)
master_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-9, "CPXPARAM_Threads" => 4)
oracle_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)
dcglp_solver_param = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9)
dcglp_param = DcglpParam(
    time_limit = 200.0,
    gap_tolerance = 1e-3,
    halt_limit = 3,
    iter_limit = 15,
    verbose = true
)
user_cb_param = UserCallbackParam(frequency=500)

# Create traditional data for MIP reference
# model = read_from_file("MIBP_domain/leo1/leo1.mps.gz") # "MIBP_domain/leo1/leo1.mps.gz"

model = read_from_file("MIBP_bound/neos-2629914-sudost/neos-2629914-sudost.mps.gz") # "MIBP_domain/neos-2629914-sudost/neos-2629914-sudost.mps.gz"

# @testset "Classical oracle" begin
#     # @testset "NoSeq" begin
#     #     data, master, typical_oracle = auto_decompose(
#     #         model;
#     #         master_solver_param = master_solver_param,
#     #         oracle_solver_param = oracle_solver_param
#     #     )
#     #     root_preprocessing = NoRootNodePreprocessing()
#     #     lazy_callback = LazyCallback(typical_oracle)
#     #     user_callback = NoUserCallback()
#     #     env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
#     #     log = solve!(env)
#     #     @test env.termination_status == Optimal()
#     #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
#     # end
#     # @testset "Seq" begin
#     #     @info "solving CFLP p$i - classical oracle - seq..."
#     #     data, master, typical_oracle = auto_decompose(
#     #         model;
#     #         master_solver_param = master_solver_param,
#     #         oracle_solver_param = oracle_solver_param
#     #     )
#     #     root_seq_type = BendersSeq
#     #     root_param = BendersSeqParam(;
#     #         time_limit = 200.0,
#     #         gap_tolerance = 1e-6,
#     #         verbose = false
#     #     )
#     #     root_preprocessing = RootNodePreprocessing(typical_oracle, root_seq_type, root_param)
#     #     lazy_callback = LazyCallback(typical_oracle)
#     #     user_callback = NoUserCallback()
#     #     env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
#     #     log = solve!(env)
#     #     @test env.termination_status == Optimal()
#     #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
#     # end
#     # @testset "SeqInOut" begin
#     #     data, master, typical_oracle = auto_decompose(
#     #         model;
#     #         master_solver_param = master_solver_param,
#     #         oracle_solver_param = oracle_solver_param
#     #     )
#     #     root_seq_type = BendersSeqInOut
#     #     root_param = BendersSeqInOutParam(
#     #         time_limit = 10.0,
#     #         gap_tolerance = 1e-6,
#     #         stabilizing_x = ones(data.dim_x),
#     #         α = 0.9,
#     #         λ = 0.1,
#     #         verbose = true
#     #     )
#     #     root_preprocessing = RootNodePreprocessing(typical_oracle, root_seq_type, root_param)
#     #     lazy_callback = LazyCallback(typical_oracle)
#     #     user_callback = NoUserCallback()
#     #     env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
#     #     log = solve!(env)
#     #     @test env.termination_status == Optimal()
#     #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
#     # end
# end

@testset "Disjunctive oracle" begin

    oracle_param = DisjunctiveOracleParam(
        norm = LpNorm(1.0),
        split_index_selection_rule = RandomFractional(),
        disjunctive_cut_append_rule = AllDisjunctiveCuts(),
        strengthened = true,
        add_benders_cuts_to_master = true,
        fraction_of_benders_cuts_to_master = 1.0,
        reuse_dcglp = true,
        adjust_t_to_fx = false
    )

    # @testset "NoSeq" begin
    #     data, master, disjunctive_oracle = auto_decompose(
    #         model,
    #         :disjunctive;
    #         master_solver_param = master_solver_param,
    #         oracle_param = oracle_param,
    #         typical_oracle_solver_param = oracle_solver_param,
    #         dcglp_solver_param = dcglp_solver_param,
    #         dcglp_param = dcglp_param
    #     )
    #     _, _, lazy_oracle = auto_decompose(
    #         model;
    #         master_solver_param = master_solver_param,
    #         oracle_solver_param = oracle_solver_param
    #     )
    #     root_preprocessing = NoRootNodePreprocessing()
    #     lazy_callback = LazyCallback(lazy_oracle)
    #     user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
    #     env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
    #     log = solve!(env)
    #     @test env.termination_status == Optimal()
    #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
    # end
    # @testset "Seq" begin
    #     @info "solving CFLP p$i - disjunctive oracle - seq..."
    #     data, master, disjunctive_oracle = auto_decompose(
    #         model,
    #         :disjunctive;
    #         master_solver_param = master_solver_param,
    #         oracle_param = oracle_param,
    #         typical_oracle_solver_param = oracle_solver_param,
    #         dcglp_solver_param = dcglp_solver_param,
    #         dcglp_param = dcglp_param
    #     )
    #     _, _, lazy_oracle = auto_decompose(
    #         model;
    #         master_solver_param = master_solver_param,
    #         oracle_solver_param = oracle_solver_param
    #     )
    #     root_seq_type = BendersSeq
    #     root_param = BendersSeqParam(;
    #         time_limit = 200.0,
    #         gap_tolerance = 1e-6,
    #         verbose = false
    #     )
    #     root_preprocessing = RootNodePreprocessing(lazy_oracle, root_seq_type, root_param)
    #     lazy_callback = LazyCallback(lazy_oracle)
    #     user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
    #     env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
    #     log = solve!(env)
    #     @test env.termination_status == Optimal()
    #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
    # end
    @testset "SeqInOut" begin
        data, master, disjunctive_oracle = auto_decompose(
            model,
            :disjunctive;
            master_solver_param = master_solver_param,
            oracle_param = oracle_param,
            typical_oracle_solver_param = oracle_solver_param,
            dcglp_solver_param = dcglp_solver_param,
            dcglp_param = dcglp_param
        )
        _, _, lazy_oracle = auto_decompose(
            model;
            master_solver_param = master_solver_param,
            oracle_solver_param = oracle_solver_param
        )
        root_seq_type = BendersSeqInOut
        root_param = BendersSeqInOutParam(
            time_limit = 10.0,
            gap_tolerance = 1e-6,
            stabilizing_x = ones(data.dim_x),
            α = 0.9,
            λ = 0.1,
            verbose = true
        )
        root_preprocessing = RootNodePreprocessing(lazy_oracle, root_seq_type, root_param)
        lazy_callback = LazyCallback(lazy_oracle)
        user_callback = UserCallback(disjunctive_oracle; params=user_cb_param)
        env = BendersBnB(data, master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
        log = solve!(env)
        @test env.termination_status == Optimal()
        @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
    end
end
