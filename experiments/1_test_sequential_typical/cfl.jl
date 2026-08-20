using BendersX
using CSV
using DataFrames
using Test
using JuMP
using CPLEX
include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))

@testset verbose = true "CFLP Sequential Benders Tests" begin
    reference_path = normpath(joinpath(@__DIR__, "..", "reference_objectives", "cflp.csv"))
    reference_df = DataFrame(CSV.File(reference_path))
    @assert nrow(reference_df) == length(unique(reference_df.instance_name)) "Duplicate CFLP reference objectives found in $(reference_path)"
    reference_objectives = Dict(String(row.instance_name) => Float64(row.objective_value) for row in eachrow(reference_df))
    instances = setdiff(1:71, [67])

    model_configs = [
        (
            name = "Classical",
            build = data -> begin
                master = Master(
                    data;
                    model = update_master_model!,
                    optimizer = mip_optimizer,
                )
                oracle = ClassicalOracle(
                    data,
                    master;
                    model = update_sub_model!,
                    optimizer = optimizer,
                )
                (master, oracle)
            end,
        ),
        (
            name = "Unified",
            build = data -> begin
                master = Master(
                    data;
                    model = update_master_model!,
                    optimizer = mip_optimizer,
                )
                oracle = UnifiedOracle(
                    data,
                    master;
                    model = update_sub_model!,
                    optimizer = optimizer,
                )
                (master, oracle)
            end,
        ),
        (
            name = "Classical + GBC",
            build = data -> begin
                master = Master(
                    data;
                    model = update_master_model!,
                    optimizer = mip_optimizer,
                )
                oracle = ClassicalOracle(
                    data,
                    master;
                    model = update_sub_gbc_model!,
                    optimizer = optimizer,
                )
                (master, oracle)
            end,
        ),
        (
            name = "Knapsack",
            build = data -> begin
                master = Master(
                    data;
                    model = update_master_model!,
                    optimizer = mip_optimizer,
                )
                oracle = CFLKnapsackOracle(
                        data, 
                        master; 
                        model = update_sub_model!, 
                        optimizer = optimizer
                )
                (master, oracle)
            end,
        ),
        (
            name = "Knapsack + GBC",
            build = data -> begin
                master = Master(
                    data;
                    model = update_master_model!,
                    optimizer = mip_optimizer,
                )
                oracle = CFLKnapsackOracle(
                        data, 
                        master; 
                        model = update_sub_gbc_model!, 
                        optimizer = optimizer
                )
                (master, oracle)
            end,
        ),
        (
            name = "Pareto",
            build = data -> begin
                master = Master(
                    data;
                    model = update_master_model!,
                    optimizer = mip_optimizer,
                )
                oracle = ParetoOracle(
                    data, 
                    master, 
                    ParetoOracleParam(fill(1.0, data.n_facilities)); 
                    model = update_sub_model!, 
                    optimizer = optimizer
                )
                (master, oracle)
            end,
        ),
    ]

    preprocessing_configs = [
        (
            name = "None",
            build = (data, oracle) -> NoPreprocessing(),
        ),
        (
            name = "LP using Seq",
            build = (data, oracle) -> LPRelaxationPreprocessing(
                oracle;
                param = BendersSeqParam(
                    time_limit = 200.0,
                    gap_tolerance = 1e-9,
                    verbose = false,
                ),
            ),
        ),
        (
            name = "LP using In-Out",
            build = (data, oracle) -> LPRelaxationPreprocessing(
                oracle;
                seq_env_type = BendersSeqInOut,
                param = BendersSeqInOutParam(
                    time_limit = 300.0,
                    gap_tolerance = 1e-9,
                    stabilizing_x = ones(data.n_facilities),
                    α = 0.9,
                    λ = 0.1,
                    verbose = false,
                ),
            ),
        ),
    ]

    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data
            data = read_cflp_benchmark_data("p$i")

            # Loop parameters
            benders_param = BendersSeqParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = false
                        )

            instance_name = "p$i"
            @assert haskey(reference_objectives, instance_name) "Missing CFLP reference objective for $(instance_name) in $(reference_path)"
            mip_opt_val = reference_objectives[instance_name]

            for model_config in model_configs
                for preprocessing_config in preprocessing_configs
                    @testset "$(model_config.name) / $(preprocessing_config.name)" begin
                        @info "solving CFLP p$i - env: BendersSeq - oracle: $(model_config.name)  - preprocessing: $(preprocessing_config.name)..."
                        master, oracle = model_config.build(data)
                        preprocessing = preprocessing_config.build(data, oracle)

                        env = BendersSeq(
                            master,
                            oracle;
                            param = benders_param,
                            preprocessing = preprocessing,
                        )

                        solve!(env)

                        @test env.termination_status == Optimal()
                        @test isapprox(mip_opt_val, env.obj_value, atol = 1e-5)
                    end
                end
            end

            # @testset "Classic oracle" begin
            #     @info "solving CFLP p$i - classical oracle - seq..."
            #     master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #     oracle = ClassicalOracle(data, master; model = update_sub_model!, optimizer = optimizer)
            #     env = BendersSeq(master, oracle; param = benders_param)
            #     log = solve!(env)
            #     @test env.termination_status == Optimal()
            #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            # end

            # @testset "Knapsack oracle" begin
            #     @info "solving CFLP p$i - knapsack oracle - seq..."
            #     master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #     oracle = CFLKnapsackOracle(data, master; model = update_sub_model!, optimizer = optimizer)
            #     env = BendersSeq(master, oracle; param = benders_param)
            #     log = solve!(env)
            #     @test env.termination_status == Optimal()
            #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            # end

            # @testset "Unified oracle" begin
            #     @info "solving CFLP p$i - unified oracle - seq..."
            #     master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #     oracle = UnifiedOracle(data, master; model = update_sub_model!, optimizer = optimizer)
            #     env = BendersSeq(master, oracle; param = benders_param)
            #     log = solve!(env)
            #     @test env.termination_status == Optimal()
            #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            # end

            # @testset "Classic oracle with GBC" begin
            #     @info "solving CFLP p$i - classical oracle with GBC - seq..."
            #     master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #     oracle = ClassicalOracle(data, master; model = update_sub_gbc_model!, optimizer = optimizer)
            #     env = BendersSeq(master, oracle; param = benders_param)
            #     log = solve!(env)
            #     @test env.termination_status == Optimal()
            #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            # end

            # @testset "Knapsack oracle with GBC" begin
            #     @info "solving CFLP p$i - knapsack oracle with GBC - seq..."
            #     master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #     oracle = CFLKnapsackOracle(data, master; model = update_sub_gbc_model!, optimizer = optimizer)
            #     env = BendersSeq(master, oracle; param = benders_param)
            #     log = solve!(env)
            #     @test env.termination_status == Optimal()
            #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            # end

            # @testset "Pareto oracle" begin
            #     @info "solving CFLP p$i - pareto oracle - seq..."
            #     master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #     param = ParetoOracleParam(fill(1.0, data.n_facilities))
            #     oracle = ParetoOracle(data, master, param; model = update_sub_model!, optimizer = optimizer)
            #     env = BendersSeq(master, oracle; param = benders_param)
            #     log = solve!(env)
            #     @test env.termination_status == Optimal()
            #     @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            # end
        end
    end
end
