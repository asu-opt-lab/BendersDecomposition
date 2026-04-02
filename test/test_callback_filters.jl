using Test
using Logging
using JuMP
using HiGHS
using BendersX

@testset "Callback metadata filters" begin
    @testset "User callback filter evaluation" begin
        result = BendersX._evaluate_user_callback_filters(UserCallbackParam())
        @test result.process_node
        @test !result.missing_node_count
        @test !result.missing_depth

        params = UserCallbackParam(node_count = 3)
        @test !BendersX._evaluate_user_callback_filters(params; node_count = 2).process_node
        @test BendersX._evaluate_user_callback_filters(params; node_count = 3).process_node
        @test BendersX._evaluate_user_callback_filters(params; node_count = 5).process_node

        params = UserCallbackParam(depth = 4)
        @test !BendersX._evaluate_user_callback_filters(params; node_depth = 3).process_node
        @test BendersX._evaluate_user_callback_filters(params; node_depth = 4).process_node
        @test BendersX._evaluate_user_callback_filters(params; node_depth = 6).process_node

        params = UserCallbackParam(node_count = 3, depth = 4)
        @test !BendersX._evaluate_user_callback_filters(params; node_count = 2, node_depth = 5).process_node
        @test !BendersX._evaluate_user_callback_filters(params; node_count = 3, node_depth = 2).process_node
        @test BendersX._evaluate_user_callback_filters(params; node_count = 3, node_depth = 4).process_node

        result = BendersX._evaluate_user_callback_filters(params; node_count = nothing, node_depth = 5)
        @test result.process_node
        @test result.missing_node_count
        @test !result.missing_depth

        result = BendersX._evaluate_user_callback_filters(params; node_count = nothing, node_depth = 2)
        @test !result.process_node
        @test result.missing_node_count
        @test !result.missing_depth

        result = BendersX._evaluate_user_callback_filters(params; node_count = 4, node_depth = nothing)
        @test result.process_node
        @test !result.missing_node_count
        @test result.missing_depth
    end

    @testset "Unsupported metadata warnings" begin
        logs, _ = Test.collect_test_logs(min_level = Logging.Warn) do
            BendersX._warn_ignored_user_callback_filters("HiGHS"; missing_node_count = true, missing_depth = true)
        end

        @test length(logs) == 1
        @test occursin("node_count", string(logs[1].message))
        @test occursin("depth", string(logs[1].message))
        @test occursin("HiGHS", string(logs[1].message))

        io = IOBuffer()
        logger = ConsoleLogger(io, Logging.Warn)
        with_logger(logger) do
            BendersX._warn_ignored_user_callback_filters("HiGHS"; missing_node_count = true)
            BendersX._warn_ignored_user_callback_filters("HiGHS"; missing_node_count = true)
        end
        output = String(take!(io))
        @test length(collect(eachmatch(r"Ignoring unsupported user callback filter\(s\)", output))) == 1
    end

    @testset "Callback metadata defaults" begin
        model = Model(HiGHS.Optimizer)
        @test BendersX.callback_node_count(nothing, model) === nothing
        @test BendersX.callback_node_depth(nothing, model) === nothing

        state = BendersX.BendersBnBState()
        @test BendersX._record_callback_node!(state, nothing) === state
        @test state.node == 0
        @test BendersX._record_callback_node!(state, 7) === state
        @test state.node == 7
    end
end
