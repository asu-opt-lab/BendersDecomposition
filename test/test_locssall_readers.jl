using Test
using BendersX

@testset "LOCSSALL benchmark readers" begin
    mktempdir() do dir
        instance_path = joinpath(dir, "pwrap")
        write(instance_path, """
        3 7
        10 100
        20 200
        30 300
        1 2 3
        4 5
        6 7
        11 12 13 14
        15 16 17
        21 22
        23 24 25 26 27
        31 32 33
        34 35 36 37
        """)

        cflp_data = read_cflp_benchmark_data("pwrap"; filepath=dir)
        @test cflp_data.n_facilities == 3
        @test cflp_data.n_customers == 7
        @test cflp_data.capacities == [10.0, 20.0, 30.0]
        @test cflp_data.fixed_costs == [100.0, 200.0, 300.0]
        @test cflp_data.demands == collect(1.0:7.0)
        @test cflp_data.costs == [
            11.0 12.0 13.0 14.0 15.0 16.0 17.0
            21.0 22.0 23.0 24.0 25.0 26.0 27.0
            31.0 32.0 33.0 34.0 35.0 36.0 37.0
        ]

        uflp_data = read_uflp_benchmark_data("pwrap"; filepath=dir)
        @test uflp_data.n_facilities == 3
        @test uflp_data.n_customers == 7
        @test uflp_data.fixed_costs == [100.0, 200.0, 300.0]
        @test uflp_data.demands == collect(1.0:7.0)
        @test uflp_data.costs == cflp_data.costs
    end
end
