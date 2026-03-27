using BendersX
import BendersX
using Test
using JuMP

@testset verbose = true "Constraint Scalarization" begin

    # ------------------------------------------------------------------ #
    # _insert_suffix                                                       #
    # ------------------------------------------------------------------ #
    @testset "_insert_suffix" begin
        # No brackets: suffix appended at end
        @test BendersX._insert_suffix("demand", "_lb") == "demand_lb"

        # Single index bracket
        @test BendersX._insert_suffix("demand[1]", "_lb") == "demand_lb[1]"

        # Multi-index bracket
        @test BendersX._insert_suffix("flow[i,j]", "_ub") == "flow_ub[i,j]"

        # Name starts with bracket (no prefix text)
        @test BendersX._insert_suffix("[1]", "_lb") == "_lb[1]"
    end

    # ------------------------------------------------------------------ #
    # MOI.Interval  →  GreaterThan + LessThan                             #
    # ------------------------------------------------------------------ #
    @testset "Interval: basic split" begin
        model = JuMP.Model()
        @variable(model, x)
        @variable(model, y)
        @constraint(model, c, 2.0 <= x + y <= 8.0)

        n_before = length(all_constraints(model, include_variable_in_set_constraints=true))
        BendersX._scalarize_constraints!(model)
        n_after = length(all_constraints(model, include_variable_in_set_constraints=true))

        @test n_after == n_before + 1

        # Original constraint removed
        @test !is_valid(model, c)

        lb_con = constraint_by_name(model, "c_lb")
        ub_con = constraint_by_name(model, "c_ub")

        # Correct constraint sense
        @test constraint_object(lb_con).set isa MOI.GreaterThan
        @test constraint_object(ub_con).set isa MOI.LessThan

        # Correct RHS values
        @test normalized_rhs(lb_con) == 2.0
        @test normalized_rhs(ub_con) == 8.0
    end

    @testset "Interval: vector of interval constraints" begin
        model = JuMP.Model()
        @variable(model, x)
        # Indexed name produces "demand[1]"
        c = @constraint(model, demand[i=1:3], i <= 2x <= i+1)

        n_before = length(all_constraints(model, include_variable_in_set_constraints=true))
        BendersX._scalarize_constraints!(model)
        n_after = length(all_constraints(model, include_variable_in_set_constraints=true))

        @test n_after == n_before + 3

        lb_con = constraint_by_name(model, "demand_lb[3]")
        ub_con = constraint_by_name(model, "demand_ub[3]")
        @test lb_con !== nothing
        @test ub_con !== nothing
        @test normalized_rhs(lb_con) == 3.0
        @test normalized_rhs(ub_con) == 4.0
    end

    @testset "Interval: unnamed constraint produces no name" begin
        model = JuMP.Model()
        @variable(model, x[1:3])
        @variable(model, y[1:5])
        @constraint(model, [i=1:3, j=i:5], 1.0 <= x[i] + y[j] <= 5.0)   # no name given

        n_before = length(all_constraints(model, include_variable_in_set_constraints=true))
        BendersX._scalarize_constraints!(model)
        n_after = length(all_constraints(model, include_variable_in_set_constraints=true))

        @test n_after == n_before + 12

        for con in all_constraints(model, AffExpr, MOI.GreaterThan{Float64})
            @test isempty(name(con))
            @test normalized_rhs(con) == 1.0
        end
        for con in all_constraints(model, AffExpr, MOI.LessThan{Float64})
            @test isempty(name(con))
            @test normalized_rhs(con) == 5.0
        end
    end

    # ------------------------------------------------------------------ #
    # MOI.Zeros  →  n × EqualTo(0)                                        #
    # (created by:  @constraint(model, A*x == b) )                        #
    # ------------------------------------------------------------------ #
    @testset "Zeros: Ax == b split" begin
        model = JuMP.Model()
        @variable(model, x[1:3])
        A = [1.0 2.0 3.0; 4.0 5.0 6.0]
        b = [5.0, 6.0]
        @constraint(model, eq_con, A * x == b)

        n_before = length(all_constraints(model, include_variable_in_set_constraints=true))
        BendersX._scalarize_constraints!(model)
        n_after = length(all_constraints(model, include_variable_in_set_constraints=true))
        @test n_after == n_before + 1

        # Original vector constraint removed
        @test !is_valid(model, eq_con)

        eq_con_1 = constraint_by_name(model, "eq_con_1")
        eq_con_2 = constraint_by_name(model, "eq_con_2")
        
        # Correct constraint sense
        @test constraint_object(eq_con_1).set isa MOI.EqualTo
        @test constraint_object(eq_con_2).set isa MOI.EqualTo

        # Correct RHS values
        @test normalized_rhs(eq_con_1) == 5.0
        @test normalized_rhs(eq_con_2) == 6.0
    end

    @testset "Zeros: Ax .== b unchanged" begin
        model = JuMP.Model()
        @variable(model, x[1:3])
        A = [1.0 3.0 5.0; 6.0 7.0 8.0]
        b = [1.0, 2.0]
        @constraint(model, eq_con, A * x .== b)

        n_before = length(all_constraints(model, include_variable_in_set_constraints=true))
        BendersX._scalarize_constraints!(model)
        n_after = length(all_constraints(model, include_variable_in_set_constraints=true))
        @test n_after == n_before

        # Original vector constraint maintained
        @test all(is_valid.(model, eq_con))
    end

    # ------------------------------------------------------------------ #
    # MOI.Nonnegatives  →  n × GreaterThan(0)                             #
    # (created by:  @constraint(model, A*x >= b) )                        #
    # ------------------------------------------------------------------ #
    @testset "Nonnegatives: named Ax >= b split" begin
        model = JuMP.Model()
        @variable(model, x[1:3])
        A = [1.0 2.0 3.0; 4.0 5.0 6.0]
        b = [2.0, 3.0]
        @constraint(model, nn_con, A * x >= b)

        BendersX._scalarize_constraints!(model)

        con1 = constraint_by_name(model, "nn_con_1")
        con2 = constraint_by_name(model, "nn_con_2")
        @test con1 !== nothing
        @test con2 !== nothing
        @test normalized_rhs(con1) == 2.0
        @test normalized_rhs(con2) == 3.0
    end

    @testset "Nonnegatives: named Ax .>= b unchanged" begin
        model = JuMP.Model()
        @variable(model, x[1:3])
        A = [1.0 2.0 3.0; 4.0 5.0 6.0]
        b = [2.0, 3.0]
        @constraint(model, nn_con, A * x .>= b)

        n_before = length(all_constraints(model, include_variable_in_set_constraints=true))
        BendersX._scalarize_constraints!(model)
        n_after = length(all_constraints(model, include_variable_in_set_constraints=true))
        @test n_after == n_before

        # Original vector constraint maintained
        @test all(is_valid.(model, nn_con))
    end

    # ------------------------------------------------------------------ #
    # MOI.Nonpositives  →  n × LessThan(0)                                #
    # (created by:  @constraint(model, A*x <= b) )                        #
    # ------------------------------------------------------------------ #
    @testset "Nonpositives: Ax <= b" begin
        model = JuMP.Model()
        @variable(model, x[1:3])
        A = [1.0 2.0 3.0; 4.0 5.0 6.0]
        b = [4.0, 5.0]
        @constraint(model, np_con, A * x <= b)

        BendersX._scalarize_constraints!(model)

        @test !is_valid(model, np_con)

        con1 = constraint_by_name(model, "np_con_1")
        con2 = constraint_by_name(model, "np_con_2")
        @test con1 !== nothing
        @test con2 !== nothing
        @test normalized_rhs(con1) == 4.0
        @test normalized_rhs(con2) == 5.0
    end

    # ------------------------------------------------------------------ #
    # No-op: already-scalar constraints are left untouched                 #
    # ------------------------------------------------------------------ #
    @testset "No-op: scalar constraints unchanged" begin
        model = JuMP.Model()
        @variable(model, x)
        @variable(model, y)
        @constraint(model, c1, x + y >= 1.0)
        @constraint(model, c2, x - y <= 5.0)
        @constraint(model, c3, x + y == 3.0)

        n_before = length(all_constraints(model, include_variable_in_set_constraints=true))
        BendersX._scalarize_constraints!(model)
        n_after = length(all_constraints(model, include_variable_in_set_constraints=true))
        @test n_before == n_after

        @test is_valid(model, c1)
        @test is_valid(model, c2)
        @test is_valid(model, c3)
        @test n_before == n_after
    end

    # ------------------------------------------------------------------ #
    # Mixed model: all four non-scalar types + pre-existing scalar         #
    # ------------------------------------------------------------------ #
    @testset "Mixed model: all types scalarized together" begin
        model = JuMP.Model()
        @variable(model, x[1:3] >= 2)
        @variable(model, sqrt(i) <= y[i = 1:3] <= i^2)
        @variable(model, z[i = 2:3, j = 1:2:3, ["red", "blue"]] >= 0)
        @variable(model, u[i = 1:2, j = i:3] >= 0)
        @variable(model, v[i = 1:9; mod(i, 3) == 0] >= 0)

        # constraints
        @constraint(model, c2, sum(i * y[i] for i in 1:3) == 3)
        @constraint(model, c3[i=1:3, j=i:3], x[i] <= j)
        @constraints(model, begin
            0.5 <= 2x[1] <= 1
            c4, x[1] >= -1
        end)
        A = [1.0 2.0 3.0; 4.0 5.0 6.0]
        b = [5.0, 6.0]
        @constraint(model, scalar_c,  x[1] + x[2] >= 0.0)
        @constraint(model, interval_c, 1.0 <= x[1] + x[2] <= 4.0)
        @constraint(model, zeros_c,   A * x == b)
        @constraint(model, nonneg_c,  A * x >= b)
        @constraint(model, nonpos_c,  A * x <= b)

        n_before = length(all_constraints(model, include_variable_in_set_constraints=true))
        BendersX._scalarize_constraints!(model)
        n_after = length(all_constraints(model, include_variable_in_set_constraints=true))
        @test n_after == n_before + 5

        # Non-scalar originals are gone
        @test !is_valid(model, interval_c)
        @test !is_valid(model, zeros_c)
        @test !is_valid(model, nonneg_c)
        @test !is_valid(model, nonpos_c)

        # Pre-existing scalar constraint survives
        @test is_valid(model, c2)
        @test all(is_valid.(model, c3))
        @test is_valid(model, c4)
        @test is_valid(model, scalar_c)

        # All remaining constraints are scalar
        for con in all_constraints(model, include_variable_in_set_constraints=true)
            set = constraint_object(con).set
            @test (set isa MOI.GreaterThan || set isa MOI.LessThan || set isa MOI.EqualTo)
        end
    end

    # ------------------------------------------------------------------ #
    # Edge case: empty model                                               #
    # ------------------------------------------------------------------ #
    @testset "Empty model: no error" begin
        model = JuMP.Model()
        @test_nowarn BendersX._scalarize_constraints!(model)
    end

end
