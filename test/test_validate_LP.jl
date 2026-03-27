using BendersX
import BendersX
import BendersX: UnsupportedModelException
using Test
using JuMP

function err_msg_test(model)
    err = try
        BendersX._validate_lp_compatibility(model)
        nothing
    catch e
        e
    end

    return err
end

@testset verbose = true "LP Compatibility Test" begin
    
    @testset "Nonlinear Objective" begin
        model = JuMP.Model()
        @variable(model, x)
        @objective(model, Min, x^2)
        
        err = err_msg_test(model)
        @test err isa UnsupportedModelException
        @test occursin("Unsupported objective function type", err.msg)
    end   

    @testset "Binary variable" begin
        model = JuMP.Model()
        @variable(model, x, Bin)
        
        err = err_msg_test(model)
        @test err isa UnsupportedModelException
        @test occursin("Discontinuous variables are not allowed.", err.msg)
    end
    
    @testset "Integer variable" begin
        model = JuMP.Model()
        @variable(model, x)
        @variable(model, y, Int)
        
        err = err_msg_test(model)
        @test err isa UnsupportedModelException
        @test occursin("Discontinuous variables are not allowed.", err.msg)
    end

    @testset "Semi-continuous variable" begin
        model = JuMP.Model()
        @variable(model, x in Semicontinuous(1.5, 3.5))
        
        err = err_msg_test(model)
        @test err isa UnsupportedModelException
        @test occursin("Discontinuous variables are not allowed.", err.msg)
    end

    @testset "PSD variable" begin
        model = JuMP.Model()
        @variable(model, x[1:2, 1:2], PSD)
        
        err = err_msg_test(model)
        @test err isa UnsupportedModelException
    end

    @testset "SOC constraint" begin
        model = JuMP.Model()
        @variable(model, t)
        @variable(model, x[1:2])
        @constraint(model, [t; x] in SecondOrderCone())
        
        err = err_msg_test(model)
        @test err isa UnsupportedModelException
    end

    @testset "SOS constraint" begin
        model = JuMP.Model()
        @variable(model, x[1:3])
        @constraint(model, x in SOS1())
        
        err = err_msg_test(model)
        @test err isa UnsupportedModelException
    end

    @testset "Quadratic constraint" begin
        model = JuMP.Model()
        @variable(model, x[i=1:2])
        @variable(model, t >= 0)
        @constraint(model, my_q, x[1]^2 + x[2]^2 <= t^2)
        
        err = err_msg_test(model)
        @test err isa UnsupportedModelException
    end

    @testset "Vectorized constraint" begin
        model = JuMP.Model()
        @variable(model, x[i=1:2])
        A = [1 2; 3 4]
        b = [5, 6]
        @constraint(model, con_vector, A * x == b)
        @constraint(model, con_scalar, A * x .== b)
        @constraint(model, A * x <= b)
        @constraint(model, A * x .<= b)
        @constraint(model, A * x >= b)
        @constraint(model, A * x .>= b)
        
        err = err_msg_test(model)
        @test err == nothing 
    end

    @testset "Allowed model" begin
        model = JuMP.Model()
    
        # various continuous variable containers
        @variable(model, x[1:3] >= 2)
        @variable(model, sqrt(i) <= y[i = 1:3] <= i^2)
        @variable(model, z[i = 2:3, j = 1:2:3, ["red", "blue"]] >= 0)
        @variable(model, u[i = 1:2, j = i:3])
        @variable(model, v[i = 1:9; mod(i, 3) == 0])
        @variable(model, w[1:2, 1:2], set = SymmetricMatrixSpace()) # symmetric variable

        # constraints
        @constraint(model, c1, 4 <= 2 * x[2] <= 5)
        @constraint(model, c2, sum(i * y[i] for i in 1:3) == 3)
        @constraint(model, c3[i=1:3, j=i:3], x[i] <= j)
        @constraints(model, begin
           2x[1] <= 1
           c4, x[1] >= -1
        end)

        err = err_msg_test(model)
        @test err == nothing
    end
end
    
