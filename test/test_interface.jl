using Test
using BendersBase
using JuMP

@testset "BendersBase copy_variables!" begin
    
    @testset "VariableRef" begin
        master_model = Model()
        @variable(master_model, u[1:10])
        @variable(master_model, v)
        x = (u = u, v = v)

        sub_model = Model()

        x_copy = copy_variables!(sub_model, x)

        function build_sub_model(model; u,v)
            @variable(model, y >= 0)
            @constraint(model, sum(u)+y == 1)
            @constraint(model, sum(v) == 1)
        end

        build_sub_model(sub_model; x_copy...)
        print(sub_model)
    end
    @testset "Array Container" begin
        @testset "1D Array" begin
            master_model = Model()
            @variable(master_model, u[1:10])
            @variable(master_model, v[1:2])
            x = (u = u, v = v)

            sub_model = Model()

            x_copy = copy_variables!(sub_model, x)

            function build_sub_model(model; u,v)
                @variable(model, y >= 0)
                @constraint(model, sum(u)+y == 1)
                @constraint(model, sum(v) == 1)
                @constraint(model, v[2]+y == 1)
            end

            build_sub_model(sub_model; x_copy...)
            print(sub_model)
        end
        @testset "2D Array" begin
            master_model = Model()
            @variable(master_model, u[1:10])
            @variable(master_model, v[1:3,1:4])
            x = (u = u, v = v)

            sub_model = Model()

            x_copy = copy_variables!(sub_model, x)

            function build_sub_model(model; u,v)
                @variable(model, y >= 0)
                @constraint(model, sum(u)+y == 1)
                @constraint(model, sum(v) == 1)
                @constraint(model, v[3,4]+y == 1)
            end

            build_sub_model(sub_model; x_copy...)
            print(sub_model)
        end
    end
    @testset "DenseAxisArray Container" begin
        @testset "1D Array" begin
            master_model = Model()
            @variable(master_model, u[1:10])
            @variable(master_model, v[2:3])
            x = (u = u, v = v)

            sub_model = Model()

            x_copy = copy_variables!(sub_model, x)

            function build_sub_model(model; u,v)
                @variable(model, y >= 0)
                @constraint(model, sum(u)+y == 1)
                @constraint(model, sum(v) == 1)
                @constraint(model, v[3]+y == 1)
            end

            build_sub_model(sub_model; x_copy...)
            print(sub_model)
        end
        @testset "1D Array with Symbol index" begin
            master_model = Model()
            @variable(master_model, u[1:10])
            @variable(master_model, v[[:A,:B]])
            x = (u = u, v = v)

            sub_model = Model()

            x_copy = copy_variables!(sub_model, x)

            function build_sub_model(model; u,v)
                @variable(model, y >= 0)
                @constraint(model, sum(u)+y == 1)
                @constraint(model, sum(v) == 1)
                @constraint(model, v[:A]+y == 1)
            end

            build_sub_model(sub_model; x_copy...)
            print(sub_model)
        end
        @testset "High dimensional array with Symbol index" begin
            master_model = Model()
            @variable(master_model, u[1:10])
            @variable(master_model, v[1:2, 2:3, [:A,:B], 4:5])
            x = (u = u, v = v)

            sub_model = Model()

            x_copy = copy_variables!(sub_model, x)

            function build_sub_model(model; u,v)
                @variable(model, y >= 0)
                @constraint(model, sum(u)+y == 1)
                @constraint(model, sum(v) == 1)
                @constraint(model, v[2,3,:A,4]+y == 1)
            end

            build_sub_model(sub_model; x_copy...)
            print(sub_model)
        end
    end
    @testset "SparseAxisArrays Container" begin
        @testset "Example 1" begin
            master_model = Model()
            @variable(master_model, u[1:10])
            @variable(master_model, v[i=1:2, j=i:2, k=j:4])
            x = (u = u, v = v)

            sub_model = Model()

            x_copy = copy_variables!(sub_model, x)

            function build_sub_model(model; u,v)
                @variable(model, y >= 0)
                @constraint(model, sum(u)+y == 1)
                @constraint(model, sum(v) == 1)
                @constraint(model, v[1,2,3]+y == 1)
            end

            build_sub_model(sub_model; x_copy...)
            print(sub_model)
        end
        @testset "Example 2" begin
            master_model = Model()
            @variable(master_model, u[1:10])
            @variable(master_model, v[i=1:4; mod(i, 2)==0])
            x = (u = u, v = v)

            sub_model = Model()

            x_copy = copy_variables!(sub_model, x)

            function build_sub_model(model; u,v)
                @variable(model, y >= 0)
                @constraint(model, sum(u)+y == 1)
                @constraint(model, sum(v) == 1)
                @constraint(model, v[2]+y == 1)
            end

            build_sub_model(sub_model; x_copy...)
            print(sub_model)
        end
    end
end