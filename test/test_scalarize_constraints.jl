using BendersX
import BendersBase
using Test
using JuMP

@testset verbose = true "Constraint Scalarization" begin

    # ------------------------------------------------------------------ #
    # _insert_suffix                                                       #
    # ------------------------------------------------------------------ #
    @testset "_insert_suffix" begin
        # No brackets: suffix appended at end
        @test BendersBase._insert_suffix("demand", "_lb") == "demand_lb"

        # Single index bracket
        @test BendersBase._insert_suffix("demand[1]", "_lb") == "demand_lb[1]"

        # Multi-index bracket
        @test BendersBase._insert_suffix("flow[i,j]", "_ub") == "flow_ub[i,j]"

        # Name starts with bracket (no prefix text)
        @test BendersBase._insert_suffix("[1]", "_lb") == "_lb[1]"
    end

    # ------------------------------------------------------------------ #
    # MOI.Interval  →  GreaterThan + LessThan                             #
    # ------------------------------------------------------------------ #
    @testset "Interval: basic split" begin
        model = JuMP.Model()
        @variable(model, x)
        @variable(model, y)
        c = @constraint(model, c, 2.0 <= x + y <= 8.0)

        _scalarize_constraints!(model)

        # Original constraint removed
        @test !is_valid(model, c)

        # Exactly one GreaterThan and one LessThan affine constraint
        ge_cons = all_constraints(model, AffExpr, MOI.GreaterThan{Float64})
        le_cons = all_constraints(model, AffExpr, MOI.LessThan{Float64})
        @test length(ge_cons) == 1
        @test length(le_cons) == 1

        # Correct RHS values
        @test normalized_rhs(first(ge_cons)) == 2.0
        @test normalized_rhs(first(le_cons)) == 8.0
    end

    @testset "Interval: named constraint suffix" begin
        model = JuMP.Model()
        @variable(model, x)
        # Indexed name produces "demand[1]"
        c = @constraint(model, demand[1], 3.0 <= 2x <= 7.0)

        _scalarize_constraints!(model)

        lb_con = constraint_by_name(model, "demand_lb[1]")
        ub_con = constraint_by_name(model, "demand_ub[1]")
        @test lb_con !== nothing
        @test ub_con !== nothing
        @test normalized_rhs(lb_con) == 3.0
        @test normalized_rhs(ub_con) == 7.0
    end

    @testset "Interval: unnamed constraint produces no name" begin
        model = JuMP.Model()
        @variable(model, x)
        @constraint(model, 1.0 <= x + 1 <= 5.0)   # no name given

        _scalarize_constraints!(model)

        for con in all_constraints(model, AffExpr, MOI.GreaterThan{Float64})
            @test isempty(name(con))
        end
        for con in all_constraints(model, AffExpr, MOI.LessThan{Float64})
            @test isempty(name(con))
        end
    end

    # ------------------------------------------------------------------ #
    # MOI.Zeros  →  n × EqualTo(0)                                        #
    # (created by:  @constraint(model, A*x == b) )                        #
    # ------------------------------------------------------------------ #
    @testset "Zeros: vector equality split" begin
        model = JuMP.Model()
        @variable(model, x[1:3])
        A = [1.0 2.0 3.0; 4.0 5.0 6.0]
        b = [5.0, 6.0]
        c = @constraint(model, eq_con, A * x == b)

        _scalarize_constraints!(model)

        # Original vector constraint removed
        @test !is_valid(model, c)

        # Two scalar EqualTo constraints
        eq_cons = all_constraints(model, AffExpr, MOI.EqualTo{Float64})
        @test length(eq_cons) == 2
    end

    @testset "Zeros: named vector equality suffix and RHS" begin
        model = JuMP.Model()
        @variable(model, x[1:3])
        A = [1.0 2.0 3.0; 4.0 5.0 6.0]
        b = [5.0, 6.0]
        @constraint(model, eq_con, A * x == b)

        _scalarize_constraints!(model)

        con1 = constraint_by_name(model, "eq_con_1")
        con2 = constraint_by_name(model, "eq_con_2")
        @test con1 !== nothing
        @test con2 !== nothing
        # Each scalar constraint encodes A[i,:]*x == b[i]
        @test normalized_rhs(con1) == 5.0
        @test normalized_rhs(con2) == 6.0
    end

    # ------------------------------------------------------------------ #
    # MOI.Nonnegatives  →  n × GreaterThan(0)                             #
    # (created by:  @constraint(model, A*x >= b) )                        #
    # ------------------------------------------------------------------ #
    @testset "Nonnegatives: vector >= split" begin
        model = JuMP.Model()
        @variable(model, x[1:3])
        A = [1.0 2.0 3.0; 4.0 5.0 6.0]
        c = @constraint(model, nn_con, A * x >= zeros(2))

        _scalarize_constraints!(model)

        @test !is_valid(model, c)
        ge_cons = all_constraints(model, AffExpr, MOI.GreaterThan{Float64})
        @test length(ge_cons) == 2
    end

    @testset "Nonnegatives: named vector >= suffix and RHS" begin
        model = JuMP.Model()
        @variable(model, x[1:3])
        A = [1.0 2.0 3.0; 4.0 5.0 6.0]
        b = [2.0, 3.0]
        @constraint(model, nn_con, A * x >= b)

        _scalarize_constraints!(model)

        con1 = constraint_by_name(model, "nn_con_1")
        con2 = constraint_by_name(model, "nn_con_2")
        @test con1 !== nothing
        @test con2 !== nothing
        @test normalized_rhs(con1) == 2.0
        @test normalized_rhs(con2) == 3.0
    end

    # ------------------------------------------------------------------ #
    # MOI.Nonpositives  →  n × LessThan(0)                                #
    # (created by:  @constraint(model, A*x <= b) )                        #
    # ------------------------------------------------------------------ #
    @testset "Nonpositives: vector <= split" begin
        model = JuMP.Model()
        @variable(model, x[1:3])
        A = [1.0 2.0 3.0; 4.0 5.0 6.0]
        c = @constraint(model, np_con, A * x <= zeros(2))

        _scalarize_constraints!(model)

        @test !is_valid(model, c)
        le_cons = all_constraints(model, AffExpr, MOI.LessThan{Float64})
        @test length(le_cons) == 2
    end

    @testset "Nonpositives: named vector <= suffix and RHS" begin
        model = JuMP.Model()
        @variable(model, x[1:3])
        A = [1.0 2.0 3.0; 4.0 5.0 6.0]
        b = [4.0, 5.0]
        @constraint(model, np_con, A * x <= b)

        _scalarize_constraints!(model)

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
        c1 = @constraint(model, c1, x + y >= 1.0)
        c2 = @constraint(model, c2, x - y <= 5.0)
        c3 = @constraint(model, c3, x + y == 3.0)

        n_before = length(all_constraints(model, include_variable_in_set_constraints=false))
        _scalarize_constraints!(model)
        n_after = length(all_constraints(model, include_variable_in_set_constraints=false))

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
        @variable(model, x[1:3])
        A = [1.0 2.0 3.0; 4.0 5.0 6.0]
        b = [5.0, 6.0]

        scalar_c  = @constraint(model, scalar_c,  x[1] + x[2] >= 0.0)
        interval_c = @constraint(model, interval_c, 1.0 <= x[1] + x[2] <= 4.0)
        zeros_c   = @constraint(model, zeros_c,   A * x == b)
        nonneg_c  = @constraint(model, nonneg_c,  A * x >= b)
        nonpos_c  = @constraint(model, nonpos_c,  A * x <= b)

        _scalarize_constraints!(model)

        # Non-scalar originals are gone
        @test !is_valid(model, interval_c)
        @test !is_valid(model, zeros_c)
        @test !is_valid(model, nonneg_c)
        @test !is_valid(model, nonpos_c)

        # Pre-existing scalar constraint survives
        @test is_valid(model, scalar_c)

        # All remaining constraints are scalar
        for con in all_constraints(model, include_variable_in_set_constraints=false)
            set = constraint_object(con).set
            @test (set isa MOI.GreaterThan || set isa MOI.LessThan || set isa MOI.EqualTo)
        end
    end

    # ------------------------------------------------------------------ #
    # Edge case: empty model                                               #
    # ------------------------------------------------------------------ #
    @testset "Empty model: no error" begin
        model = JuMP.Model()
        @test_nowarn _scalarize_constraints!(model)
    end

end
