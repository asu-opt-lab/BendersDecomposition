using Test
using Logging
using JuMP
using HiGHS
using BendersX

@testset "Callback metadata thresholds" begin
    @testset "Unsupported metadata warnings" begin
        model = Model(HiGHS.Optimizer)

        logs, _ = Test.collect_test_logs(min_level = Logging.Warn) do
            BendersX.callback_node_count(nothing, model)
            BendersX.callback_node_depth(nothing, model)
        end

        @test length(logs) == 2
        @test occursin("node_count", string(logs[1].message))
        @test occursin("HiGHS", string(logs[1].message))
        @test occursin("depth", string(logs[2].message))
        @test occursin("HiGHS", string(logs[2].message))
    end

    @testset "Unsupported metadata warnings are throttled" begin
        model = Model(HiGHS.Optimizer)
        io = IOBuffer()
        logger = ConsoleLogger(io, Logging.Warn)
        with_logger(logger) do
            BendersX.callback_node_count(nothing, model)
            BendersX.callback_node_count(nothing, model)
        end
        output = String(take!(io))
        @test length(collect(eachmatch(r"callback node count metadata", output))) == 1
    end

    @testset "Callback metadata defaults" begin
        model = Model(HiGHS.Optimizer)
        logs, _ = Test.collect_test_logs(min_level = Logging.Warn) do
            @test BendersX.callback_node_count(nothing, model) === nothing
            @test BendersX.callback_node_depth(nothing, model) === nothing
        end
        @test length(logs) == 2

        state = BendersX.BendersBnBState()
        @test BendersX.record_node_count!(state, nothing) === state
        @test state.node == 0
        @test BendersX.record_node_count!(state, 7) === state
        @test state.node == 7
    end
end
