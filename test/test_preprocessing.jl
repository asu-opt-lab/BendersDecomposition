using Test
using JuMP

struct PreprocessingTestOracle <: BendersX.AbstractOracle end

struct PreprocessingTestMaster <: BendersX.AbstractMaster
    model::Model
end

struct PreprocessingTestSplitOracle <: BendersX.AbstractSplitOracle
    param::Any
end

PreprocessingTestSplitOracle() = PreprocessingTestSplitOracle((
    split_index_selection_rule = LargestFractional(),
    disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(),
))

struct PreprocessingProbeException <: Exception end
struct ThrowingPreprocessing <: BendersX.AbstractBendersPreprocessing end

BendersX.preprocess!(::BendersX.AbstractMaster, ::ThrowingPreprocessing) =
    throw(PreprocessingProbeException())

@testset "Sequential preprocessing" begin
    @testset "LPRelaxationPreprocessing keyword constructor" begin
        oracle = PreprocessingTestOracle()
        param = BendersSeqParam(verbose = false)

        keyword = LPRelaxationPreprocessing(
            oracle;
            seq_env_type = BendersSeq,
            param = param,
        )

        @test keyword.oracle === oracle
        @test keyword.seq_env_type === BendersSeq
        @test keyword.param === param
        @test_throws MethodError LPRelaxationPreprocessing(oracle, BendersSeq, param)
    end

    @testset "BendersSeqInOut preprocessing" begin
        preprocessing = ThrowingPreprocessing()
        env = BendersSeqInOut(
            PreprocessingTestMaster(Model()),
            PreprocessingTestOracle();
            param = BendersSeqInOutParam(
                stabilizing_x = Float64[],
                verbose = false,
            ),
            preprocessing = preprocessing,
        )

        @test env.preprocessing === preprocessing
        @test_throws PreprocessingProbeException solve!(env)
    end

    @testset "SpecializedBendersSeq preprocessing" begin
        preprocessing = ThrowingPreprocessing()
        env = SpecializedBendersSeq(
            PreprocessingTestMaster(Model()),
            PreprocessingTestSplitOracle();
            param = SpecializedBendersSeqParam(verbose = false),
            preprocessing = preprocessing,
        )

        @test env.preprocessing === preprocessing
        @test_throws PreprocessingProbeException solve!(env)
    end
end
