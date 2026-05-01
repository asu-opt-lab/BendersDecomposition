using Statistics

include(joinpath(@__DIR__, "..", "src", "read_flcap.jl"))

@testset "FLCAP Reader" begin
    data = read_flcap_data("cap61")

    @test data.n_facilities == 16
    @test data.n_customers == 50
    @test data.demands[1] == 146.0
    @test isapprox(data.costs[1, 1], 46.1625; atol = 1e-9)
    @test isapprox(data.costs[1, 2], 36.8375; atol = 1e-9)

    cost_per_unit = vec(data.costs)
    @test median(cost_per_unit) < 200.0
end
