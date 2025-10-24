using Test
using BendersBase

@testset "BendersBase Types" begin
    

    @testset "Abstract Types" begin
        @test isabstracttype(AbstractBendersDecomposition)
        @test isabstracttype(AbstractMip)
        @test isabstracttype(AbstractMaster)
        @test isabstracttype(AbstractOracle)
        @test isabstracttype(AbstractBendersSeq)
        @test isabstracttype(AbstractBendersCallback)
    end
    
end