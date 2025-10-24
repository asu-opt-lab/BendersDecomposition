using Test
using BendersBase

@testset "BendersBase Types" begin
    

    @testset "Abstract Types" begin
        @test isabstracttype(AbstractNorm)
    end
    
end