using Test

struct PreprocessingTestOracle <: BendersX.AbstractOracle end

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
end
