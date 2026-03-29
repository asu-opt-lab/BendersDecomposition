using JuMP
using Test

function _capture_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

@testset "Optional Solver Extensions" begin
    cplex_err = _capture_error(() -> BendersX.cplex_optimizer())
    @test cplex_err isa ArgumentError
    @test occursin("Install and load CPLEX", sprint(showerror, cplex_err))

    gurobi_err = _capture_error(() -> BendersX.gurobi_optimizer())
    @test gurobi_err isa ArgumentError
    @test occursin("Install and load Gurobi", sprint(showerror, gurobi_err))

    missing_solver_err =
        _capture_error(() -> BendersX.assign_attributes!(Model(), Dict{String,Any}()))
    @test missing_solver_err isa ArgumentError
    @test occursin("\"solver\"", sprint(showerror, missing_solver_err))

    unsupported_solver_err =
        _capture_error(() -> BendersX.assign_attributes!(Model(), Dict("solver" => "HiGHS")))
    @test unsupported_solver_err isa ArgumentError
    @test occursin("Unsupported solver", sprint(showerror, unsupported_solver_err))

    delegated_cplex_err =
        _capture_error(() -> BendersX.assign_attributes!(Model(), Dict("solver" => "CPLEX")))
    @test delegated_cplex_err isa ArgumentError
    @test occursin("Install and load CPLEX", sprint(showerror, delegated_cplex_err))

    delegated_gurobi_err =
        _capture_error(() -> BendersX.assign_attributes!(Model(), Dict("solver" => "gurobi")))
    @test delegated_gurobi_err isa ArgumentError
    @test occursin("Install and load Gurobi", sprint(showerror, delegated_gurobi_err))
end
