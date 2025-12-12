using Test
using BendersBase

@testset "BendersBase Types" begin
    

    @testset "Abstract Types" begin
        @test isabstracttype(AbstractBendersEnv)
        @test isabstracttype(AbstractMaster)
        @test isabstracttype(AbstractOracle)
        @test isabstracttype(AbstractBendersSeq)
        @test isabstracttype(AbstractBendersBnB)
    end
    
end