using BendersX
using CSV
using DataFrames
using Test
using JuMP
using CPLEX
include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))

@testset verbose = true "UFLP Sequential Benders Tests" begin
    reference_path = normpath(joinpath(@__DIR__, "..", "reference_objectives", "uflp.csv"))
    reference_df = DataFrame(CSV.File(reference_path))
    @assert nrow(reference_df) == length(unique(reference_df.instance_name)) "Duplicate UFLP reference objectives found in $(reference_path)"
    reference_objectives = Dict(String(row.instance_name) => Float64(row.objective_value) for row in eachrow(reference_df))
    
    instances = setdiff(1:71, [67])

    # GBC-enabled subproblem customization (y[i,j] <= x[i] via GBC)
    function update_sub_gbc_model!(model::Model, data::UFLPData, scen_idx::Int; x)
        optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, MOI.Silent() => true)
        set_optimizer(model, optimizer)

        I, J = data.n_facilities, data.n_customers
        @variable(model, y[1:I, 1:J] >= 0)

        cost_demands = data.costs .* data.demands'
        @objective(model, Min, sum(cost_demands .* y))

        @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)

        # Return GBC tuple: y[i,j] <= x[i] for each j
        gbc_lhs = vec(y)
        gbc_rhs = [x[i] for j in 1:J for i in 1:I]  # j outer, i inner to match vec(y)
        gbc_sense = fill(UpperBound, I*J)
        return gbc_lhs, gbc_rhs, gbc_sense
    end

    function update_knapsack_master_model!(model::Model, data::UFLPData)
        optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)
        set_optimizer(model, optimizer)
        @variable(model, x[1:data.n_facilities], Bin)
        @variable(model, t[1:data.n_customers] >= -1e6)
        @constraint(model, sum(x) >= 2)
        @objective(model, Min, data.fixed_costs'* x + sum(t))
        return (x = x, ), t
    end

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
                    model = update_knapsack_master_model!,
                    optimizer = mip_optimizer,
                )
                oracle = UFLKnapsackOracle(
                    data;
                    param = UFLKnapsackOracleParam(
                        add_only_violated_cuts = true,
                    ),
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
                    verbose = true,
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
                    verbose = true,
                ),
            ),
        ),
    ]

    for i in instances
        @testset "Instance: p$i" begin

            # Load problem data
            data = read_uflp_benchmark_data("p$(i)")

            # Loop parameters
            benders_param = BendersSeqParam(;
                            time_limit = 200.0,
                            gap_tolerance = 1e-6,
                            verbose = true
                            )

            instance_name = "p$i"
            @assert haskey(reference_objectives, instance_name) "Missing UFLP reference objective for $(instance_name) in $(reference_path)"
            mip_opt_val = reference_objectives[instance_name]

            for model_config in model_configs
                for preprocessing_config in preprocessing_configs
                    @testset "$(model_config.name) / $(preprocessing_config.name)" begin
                        @info "solving UFLP p$i - oracle: $(model_config.name)  - preprocessing: $(preprocessing_config.name)..."
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
            #     @testset "No preprocessing" begin
            #         @info "solving UFLP p$i - classical oracle - no preprocessing..."

            #         master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #         oracle = ClassicalOracle(data, master; model = update_sub_model!, optimizer = optimizer)
            #         env = BendersSeq(master, oracle; param = benders_param)
            #         log = solve!(env)
            #         @test env.termination_status == Optimal()
            #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            #     end

            #     @testset "LP preprocessing using BendersSeq" begin
            #         @info "solving UFLP p$i - classical oracle - LP preprocessing w/ BendersSeq..."

            #         master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #         oracle = ClassicalOracle(data, master; model = update_sub_model!, optimizer = optimizer)

            #         preprocessing_seq_param = BendersSeqParam(;
            #                     time_limit = 200.0,
            #                     gap_tolerance = 1e-9,
            #                     verbose = true
            #                 )
            #         preprocessing = LPRelaxationPreprocessing(oracle; param = preprocessing_seq_param)

            #         env = BendersSeq(master, oracle; param = benders_param, preprocessing = preprocessing)
            #         log = solve!(env)
            #         @test env.termination_status == Optimal()
            #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            #     end

            #     @testset "LP preprocessing using BendersSeqInOut" begin
            #         @info "solving UFLP p$i - classical oracle - LP preprocessing w/ BendersSeqInOut..."
            #         master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #         oracle = ClassicalOracle(data, master; model = update_sub_model!, optimizer = optimizer)

            #         preprocessing_seq_type = BendersSeqInOut
            #         preprocessing_seq_param = BendersSeqInOutParam(
            #                     time_limit = 300.0,
            #                     gap_tolerance = 1e-9,
            #                     stabilizing_x = ones(data.n_facilities),
            #                     α = 0.9,
            #                     λ = 0.1,
            #                     verbose = true
            #                 )

            #         preprocessing = LPRelaxationPreprocessing(oracle; seq_env_type =  preprocessing_seq_type, param = preprocessing_seq_param)
                    
            #         env = BendersSeq(master, oracle; param = benders_param, preprocessing = preprocessing)
            #         log = solve!(env)
            #         @test env.termination_status == Optimal()
            #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            #     end
            # end

            # @testset "Unified oracle" begin
            #     @testset "No preprocessing" begin
            #         @info "solving UFLP p$i - unified oracle - no preprocessing..."
            #         master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #         oracle = UnifiedOracle(data, master; model = update_sub_model!, optimizer = optimizer)
            #         env = BendersSeq(master, oracle; param = benders_param)
            #         log = solve!(env)
            #         @test env.termination_status == Optimal()
            #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            #     end

            #     @testset "LP preprocessing using BendersSeq" begin
            #         @info "solving UFLP p$i - unified oracle - LP preprocessing w/ BendersSeq..."

            #         master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #         oracle = UnifiedOracle(data, master; model = update_sub_model!, optimizer = optimizer)

            #         preprocessing_seq_param = BendersSeqParam(;
            #                     time_limit = 200.0,
            #                     gap_tolerance = 1e-9,
            #                     verbose = true
            #                 )
            #         preprocessing = LPRelaxationPreprocessing(oracle; param = preprocessing_seq_param)

            #         env = BendersSeq(master, oracle; param = benders_param, preprocessing = preprocessing)
            #         log = solve!(env)
            #         @test env.termination_status == Optimal()
            #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            #     end

            #     @testset "LP preprocessing using BendersSeqInOut" begin
            #         @info "solving UFLP p$i - unified oracle - LP preprocessing w/ BendersSeqInOut..."
            #         master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #         oracle = UnifiedOracle(data, master; model = update_sub_model!, optimizer = optimizer)

            #         preprocessing_seq_type = BendersSeqInOut
            #         preprocessing_seq_param = BendersSeqInOutParam(
            #                     time_limit = 300.0,
            #                     gap_tolerance = 1e-9,
            #                     stabilizing_x = ones(data.n_facilities),
            #                     α = 0.9,
            #                     λ = 0.1,
            #                     verbose = true
            #                 )

            #         preprocessing = LPRelaxationPreprocessing(oracle; seq_env_type = preprocessing_seq_type, param = preprocessing_seq_param)
                    
            #         env = BendersSeq(master, oracle; param = benders_param, preprocessing = preprocessing)
            #         log = solve!(env)
            #         @test env.termination_status == Optimal()
            #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            #     end
            # end

            # @testset "Classic oracle with GBC" begin
            #     @testset "No preprocessing" begin
            #         @info "solving UFLP p$i - classical oracle with GBC - no preprocessing..."
            #         master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #         oracle = ClassicalOracle(data, master; model = update_sub_gbc_model!, optimizer = optimizer)
            #         env = BendersSeq(master, oracle; param = benders_param)
            #         log = solve!(env)
            #         @test env.termination_status == Optimal()
            #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            #     end

            #     @testset "LP preprocessing using BendersSeq" begin
            #         @info "solving UFLP p$i - classical oracle with GBC - LP preprocessing w/ BendersSeq..."

            #         master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #         oracle = ClassicalOracle(data, master; model = update_sub_gbc_model!, optimizer = optimizer)

            #         preprocessing_seq_param = BendersSeqParam(;
            #                     time_limit = 200.0,
            #                     gap_tolerance = 1e-9,
            #                     verbose = true
            #                 )
            #         preprocessing = LPRelaxationPreprocessing(oracle; param = preprocessing_seq_param)

            #         env = BendersSeq(master, oracle; param = benders_param, preprocessing = preprocessing)
            #         log = solve!(env)
            #         @test env.termination_status == Optimal()
            #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            #     end

            #     @testset "LP preprocessing using BendersSeqInOut" begin
            #         @info "solving UFLP p$i - classical oracle with GBC - LP preprocessing w/ BendersSeqInOut..."
            #         master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #         oracle = ClassicalOracle(data, master; model = update_sub_gbc_model!, optimizer = optimizer)

            #         preprocessing_seq_type = BendersSeqInOut
            #         preprocessing_seq_param = BendersSeqInOutParam(
            #                     time_limit = 300.0,
            #                     gap_tolerance = 1e-9,
            #                     stabilizing_x = ones(data.n_facilities),
            #                     α = 0.9,
            #                     λ = 0.1,
            #                     verbose = true
            #                 )

            #         preprocessing = LPRelaxationPreprocessing(oracle; seq_env_type = preprocessing_seq_type, param = preprocessing_seq_param)
                    
            #         env = BendersSeq(master, oracle; param = benders_param, preprocessing = preprocessing)
            #         log = solve!(env)
            #         @test env.termination_status == Optimal()
            #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            #     end
            # end

            # @testset "Knapsack oracle" begin

            #     function update_master_model!(model::Model, data::UFLPData)
            #         optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)
            #         set_optimizer(model, optimizer)
            #         @variable(model, x[1:data.n_facilities], Bin)
            #         @variable(model, t[1:data.n_customers] >= -1e6)
            #         @constraint(model, sum(x) >= 2)
            #         @objective(model, Min, data.fixed_costs'* x + sum(t))
            #         return (x = x, ), t
            #     end

            #     @testset "Fat version" begin

            #         @info "solving UFLP p$i - fat knapsack oracle - seq..."

            #         master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
            #         oracle = UFLKnapsackOracle(
            #             data;
            #             param = UFLKnapsackOracleParam(add_only_violated_cuts = true),
            #         )

            #         env = BendersSeq(master, oracle; param = benders_param)
            #         log = solve!(env)
            #         @test env.termination_status == Optimal()
            #         @test isapprox(mip_opt_val, env.obj_value, atol=1e-5)
            #     end

            #     # To test the slim version, pass UFLKnapsackOracleParam(slim = true) to the oracle constructor.
            # end

            # @testset "Pareto oracle" begin
            #     @info "solving UFLP p$i - pareto oracle - seq..."
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
