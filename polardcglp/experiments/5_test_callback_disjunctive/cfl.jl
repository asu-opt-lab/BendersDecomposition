using BendersX
using CSV
using DataFrames
using Test
using JuMP
using CPLEX

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))

@testset verbose = true "CFLP Disjunctive Callback Benders Tests" begin
    root_experiments_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "experiments"))
    reference_path = normpath(joinpath(root_experiments_dir, "reference_objectives", "cflp.csv"))
    reference_df = DataFrame(CSV.File(reference_path))
    @assert nrow(reference_df) == length(unique(reference_df.instance_name)) "Duplicate CFLP reference objectives found in $(reference_path)"
    reference_objectives = Dict(String(row.instance_name) => Float64(row.objective_value) for row in eachrow(reference_df))
    instances = setdiff(1:71, [67])

    for i in instances
        @testset "Instance: p$i" begin
            data = read_cflp_benchmark_data("p$i")

            benders_param = BendersBnBParam(
                time_limit = 200.0,
                gap_tolerance = 1e-6,
                verbose = false,
            )
            dcglp_optimizer = optimizer_with_attributes(
                CPLEX.Optimizer,
                "CPXPARAM_Threads" => 7,
                "CPX_PARAM_EPRHS" => 1e-9,
                "CPX_PARAM_NUMERICALEMPHASIS" => 1,
                "CPX_PARAM_EPOPT" => 1e-9,
                MOI.Silent() => true,
            )
            dcglp_param = DcglpParam(
                dcglp_optimizer;
                time_limit = 200.0,
                gap_tolerance = 1e-3,
                halt_limit = 3,
                iter_limit = 15,
                verbose = false,
            )
            user_cb_param = UserCallbackParam(frequency = 1)

            instance_name = "p$i"
            @assert haskey(reference_objectives, instance_name) "Missing CFLP reference objective for $(instance_name) in $(reference_path)"
            mip_opt_val = reference_objectives[instance_name]

            @testset "Classical oracle" begin
                oracle_param = SimplexNormDCGLPParam(
                    dcglp_param;
                    split_index_selection_rule = RandomFractional(),
                    disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                    add_benders_cuts_to_master = true,
                    fraction_of_benders_cuts_to_master = 0.5,
                    reuse_dcglp = true,

                    zero_tol = 1e-9,
                )

                @testset "NoSeq" begin
                    @info "solving SimplexNormDCGLP CFLP p$i - callback - classical/no seq"
                    master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)
                    lazy_oracle = ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer)
                    typical_oracles = [
                        ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                        ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                    ]
                    disjunctive_oracle = SimplexNormDCGLPOracle(master, typical_oracles, oracle_param)

                    root_preprocessing = NoRootNodePreprocessing()
                    lazy_callback = LazyCallback(lazy_oracle)
                    user_callback = UserCallback(disjunctive_oracle; params = user_cb_param)

                    env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                end

                @testset "Seq" begin
                    @info "solving SimplexNormDCGLP CFLP p$i - callback - classical/seq"
                    master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)
                    lazy_oracle = ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer)
                    typical_oracles = [
                        ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                        ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                    ]
                    disjunctive_oracle = SimplexNormDCGLPOracle(master, typical_oracles, oracle_param)

                    root_preprocessing = RootNodePreprocessing(
                        lazy_oracle,
                        BendersSeq,
                        BendersSeqParam(time_limit = 200.0, gap_tolerance = 1e-9, verbose = false),
                    )
                    lazy_callback = LazyCallback(lazy_oracle)
                    user_callback = UserCallback(disjunctive_oracle; params = user_cb_param)

                    env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                end

                @testset "SeqInOut" begin
                    @info "solving SimplexNormDCGLP CFLP p$i - callback - classical/seqinout"
                    master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)
                    lazy_oracle = ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer)
                    typical_oracles = [
                        ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                        ClassicalOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                    ]
                    disjunctive_oracle = SimplexNormDCGLPOracle(master, typical_oracles, oracle_param)

                    root_preprocessing = RootNodePreprocessing(
                        lazy_oracle,
                        BendersSeqInOut,
                        BendersSeqInOutParam(
                            time_limit = 300.0,
                            gap_tolerance = 1e-9,
                            stabilizing_x = ones(data.n_facilities),
                            α = 0.9,
                            λ = 0.1,
                            verbose = false,
                        ),
                    )
                    lazy_callback = LazyCallback(lazy_oracle)
                    user_callback = UserCallback(disjunctive_oracle; params = user_cb_param)

                    env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                end
            end

            @testset "Knapsack oracle" begin
                oracle_param = SimplexNormDCGLPParam(
                    dcglp_param;
                    split_index_selection_rule = RandomFractional(),
                    disjunctive_cut_append_rule = AllDisjunctiveCuts(),
                    add_benders_cuts_to_master = true,
                    fraction_of_benders_cuts_to_master = 0.5,
                    reuse_dcglp = true,

                    zero_tol = 1e-9,
                )

                @testset "NoSeq" begin
                    @info "solving SimplexNormDCGLP CFLP p$i - callback - knapsack/no seq"
                    master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)
                    lazy_oracle = CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer)
                    typical_oracles = [
                        CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                        CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                    ]
                    disjunctive_oracle = SimplexNormDCGLPOracle(master, typical_oracles, oracle_param)

                    root_preprocessing = NoRootNodePreprocessing()
                    lazy_callback = LazyCallback(lazy_oracle)
                    user_callback = UserCallback(disjunctive_oracle; params = user_cb_param)

                    env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                end

                @testset "Seq" begin
                    @info "solving SimplexNormDCGLP CFLP p$i - callback - knapsack/seq"
                    master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)
                    lazy_oracle = CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer)
                    typical_oracles = [
                        CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                        CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                    ]
                    disjunctive_oracle = SimplexNormDCGLPOracle(master, typical_oracles, oracle_param)

                    root_preprocessing = RootNodePreprocessing(
                        lazy_oracle,
                        BendersSeq,
                        BendersSeqParam(time_limit = 200.0, gap_tolerance = 1e-9, verbose = false),
                    )
                    lazy_callback = LazyCallback(lazy_oracle)
                    user_callback = UserCallback(disjunctive_oracle; params = user_cb_param)

                    env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                end

                @testset "SeqInOut" begin
                    @info "solving SimplexNormDCGLP CFLP p$i - callback - knapsack/seqinout"
                    master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)
                    lazy_oracle = CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer)
                    typical_oracles = [
                        CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                        CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
                    ]
                    disjunctive_oracle = SimplexNormDCGLPOracle(master, typical_oracles, oracle_param)

                    root_preprocessing = RootNodePreprocessing(
                        lazy_oracle,
                        BendersSeqInOut,
                        BendersSeqInOutParam(
                            time_limit = 300.0,
                            gap_tolerance = 1e-9,
                            stabilizing_x = ones(data.n_facilities),
                            α = 0.9,
                            λ = 0.1,
                            verbose = false,
                        ),
                    )
                    lazy_callback = LazyCallback(lazy_oracle)
                    user_callback = UserCallback(disjunctive_oracle; params = user_cb_param)

                    env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)
                    solve!(env)
                    @test env.termination_status == Optimal()
                    @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                end
            end
        end
    end
end
