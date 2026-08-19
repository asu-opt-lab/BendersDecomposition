struct CallbackTestDisjunctiveOracle <: BendersX.AbstractDisjunctiveOracle end

@testset "Callback configuration" begin
    oracle = CallbackTestDisjunctiveOracle()
    callback = LazyCallback(oracle)

    @test callback.oracle === oracle
end
