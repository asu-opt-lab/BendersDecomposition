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

@testset "BendersBase customize model functions" begin
    struct EmptyData <: AbstractData end
    problem = EmptyData()

    @testset "master variable container Vector{VariableRef}" begin
        function customize_master_model!(model::Model, problem::EmptyData)

            @variable(model, u[1:10], Bin)
            @variable(model, t >= -1e6)
            @constraint(model, sum(u) >= 2)
            @objective(model, Min, 1.0 * t)
            
            return (u = u, ), t
        end
        
        function customize_sub_model!(model::Model, problem::EmptyData, scen_idx::Int; u) 
        
            @variable(model, y[1:10] >= 0)
            @objective(model, Min, sum(y))
            @constraint(model, y .<= u)
        end

        master = Master(problem; customize = customize_master_model!)
        oracle = ClassicalOracle(problem, master; customize = customize_sub_model!)

        print(master.model)
        print(oracle.model)
    end
    @testset "master variable container Matrix{VariableRef}" begin
        function customize_master_model!(model::Model, problem::EmptyData)

            @variable(model, u[1:10,1:3], Bin)
            @variable(model, t >= -1e6)
            @constraint(model, sum(u) >= 2)
            @objective(model, Min, 1.0 * t)
            
            return (u = u, ), t
        end
        
        function customize_sub_model!(model::Model, problem::EmptyData, scen_idx::Int; u) 
        
            @variable(model, y[1:10] >= 0)
            @objective(model, Min, sum(y))
            @constraint(model, sum(u) == 1)
            @constraint(model, y .<= u[:,1])
        end

        master = Master(problem; customize = customize_master_model!)
        oracle = ClassicalOracle(problem, master; customize = customize_sub_model!)

        print(master.model)
        print(oracle.model)
    end
    @testset "Multiple master variable containers of varying types" begin
        function customize_master_model!(model::Model, problem::EmptyData)

            @variable(model, u[1:10], Bin) # Array
            @variable(model, v[3:10, [:A, :B]], Bin) # DenseAxisArray
            @variable(model, w[i=1:3,j=i:10], Bin) # SparseAxisArray
            @variable(model, t >= -1e6)
            @constraint(model, sum(u) >= 2)
            @objective(model, Min, 1.0 * t)
            
            return (u = u, v = v, w = w), t
        end
        
        function customize_sub_model!(model::Model, problem::EmptyData, scen_idx::Int; u, v, w) 
        
            @variable(model, y[1:10] >= 0)
            @objective(model, Min, sum(y))
            @constraint(model, y .<= u)
            @constraint(model, sum(v) == 1)
            @constraint(model, v[10,:A] == 1)
            @constraint(model, w[3,10] <= 1)
        end

        master = Master(problem; customize = customize_master_model!)
        oracle = ClassicalOracle(problem, master; customize = customize_sub_model!)

        print(master.model)
        print(oracle.model)
    end
    
end