using Test
using BendersDecomposition
using Statistics

@testset "Stochastic Capacited Facility Location Data Generation" begin
    # Test case parameters
    n_facilities = 3
    n_customers = 4
    n_scenarios = 2
    ratio = 2

    # Generate test data
    data = generate_stochastic_capacited_facility_location(
        n_facilities,
        n_customers,
        n_scenarios,
        ratio
    )

    # Basic structure tests
    @test data.n_facilities == n_facilities
    @test data.n_customers == n_customers
    @test data.n_scenarios == n_scenarios
    
    # Dimension tests
    @test length(data.capacities) == n_facilities
    @test length(data.fixed_costs) == n_facilities
    @test length(data.demands) == n_scenarios
    @test all(length.(data.demands) .== n_customers)
    @test size(data.costs) == (n_customers, n_facilities)

    # Value range tests
    @test all(data.capacities .>= 0)
    @test all(data.fixed_costs .>= 0)
    @test all(vcat(data.demands...) .>= 1)  # 确保所有需求至少为1
    @test all(data.costs .>= 0)

    # Capacity ratio test
    max_total_demand = maximum([sum(d) for d in data.demands])
    total_capacity = sum(data.capacities)
    @test total_capacity ≈ ratio * max_total_demand rtol=0.1

    # Demand variation tests
    if n_scenarios > 1
        # 检查不同场景的需求是否有变化
        demands_matrix = hcat(data.demands...)  # 将demands转换为矩阵以便计算
        @test std(sum.(data.demands)) > 0  # 确保不同场景的总需求有变化
    end
end
