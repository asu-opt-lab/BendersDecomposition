"""
UnifiedOracle Test Suite

Tests the UnifiedOracle implementation for Benders decomposition.
Compares results between ClassicalOracle and UnifiedOracle to ensure consistency.

Uses a simple assignment problem structure for testing.
"""

using BendersX
using Test
using JuMP
using HiGHS

# =============================================================================
# Test Problem: Simple Assignment Problem
# =============================================================================

struct SimpleUnifiedTestData <: AbstractData
    n_facilities::Int
    n_customers::Int
    costs::Matrix{Float64}
end

function create_unified_test_data(n_facilities::Int=3, n_customers::Int=4)
    # Use fixed seed for reproducibility
    costs = [1.0 2.0 3.0 4.0;
             5.0 1.0 2.0 3.0;
             4.0 3.0 1.0 2.0]
    return SimpleUnifiedTestData(n_facilities, n_customers, costs)
end

function customize_master_unified!(model::Model, data::SimpleUnifiedTestData)
    optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
    set_optimizer(model, optimizer)
    
    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t >= -1e6)
    @objective(model, Min, t)
    
    # At least one facility must be open
    @constraint(model, sum(x) >= 1)
    
    return (x = x, ), t
end

function solve_mip_unified_reference(data::SimpleUnifiedTestData)
    model = Model()
    optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
    set_optimizer(model, optimizer)
    
    I, J = data.n_facilities, data.n_customers
    @variable(model, x[1:I], Bin)
    @variable(model, y[1:I, 1:J] >= 0)
    
    @objective(model, Min, sum(data.costs .* y))
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= x[i])
    @constraint(model, sum(x) >= 1)
    
    optimize!(model)
    return objective_value(model)
end

# =============================================================================
# Test Scenarios
# =============================================================================

@testset verbose = true "UnifiedOracle Test Suite" begin
    
    @testset "Basic UnifiedOracle Construction" begin
        # Test that UnifiedOracle can be constructed
        data = create_unified_test_data()
        
        function customize_sub_basic!(model::Model, data::SimpleUnifiedTestData, scen_idx::Int; x)
            optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            I, J = data.n_facilities, data.n_customers
            @variable(model, y[1:I, 1:J] >= 0)
            
            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= x[i])
            
            return nothing
        end
        
        master = Master(data; customize = customize_master_unified!)
        oracle = UnifiedOracle(data, master; customize = customize_sub_basic!)
        
        # Verify oracle has all required fields
        @test !isempty(oracle.fixing_lb_constraints)
        @test !isempty(oracle.fixing_ub_constraints)
        @test oracle.objective_constraint isa ConstraintRef
        @test oracle.model isa Model
    end
    
    @testset "UnifiedOracle vs ClassicalOracle - Basic Problem" begin
        # Compare results between ClassicalOracle and UnifiedOracle
        data = create_unified_test_data()
        mip_opt = solve_mip_unified_reference(data)
        
        function customize_sub_compare!(model::Model, data::SimpleUnifiedTestData, scen_idx::Int; x)
            optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            I, J = data.n_facilities, data.n_customers
            @variable(model, y[1:I, 1:J] >= 0)
            
            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= x[i])
            
            return nothing
        end
        
        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        
        # Test ClassicalOracle
        master_classical = Master(data; customize = customize_master_unified!)
        oracle_classical = ClassicalOracle(data, master_classical; customize = customize_sub_compare!)
        env_classical = BendersSeq(master_classical, oracle_classical; param = benders_param)
        solve!(env_classical)
        
        @test env_classical.termination_status == Optimal()
        @test isapprox(mip_opt, env_classical.obj_value, atol=1e-4)
        
        # Test UnifiedOracle
        master_unified = Master(data; customize = customize_master_unified!)
        oracle_unified = UnifiedOracle(data, master_unified; customize = customize_sub_compare!)
        env_unified = BendersSeq(master_unified, oracle_unified; param = benders_param)
        solve!(env_unified)
        
        @test env_unified.termination_status == Optimal()
        @test isapprox(mip_opt, env_unified.obj_value, atol=1e-4)
        
        # Both should find the same solution
        @test isapprox(env_classical.obj_value, env_unified.obj_value, atol=1e-4)
    end
    
    @testset "UnifiedOracle Parameters" begin
        # Test that UnifiedOracleParam can be customized
        param = UnifiedOracleParam(rtol = 1e-8, atol = 1e-6, zero_tol = 1e-5)
        
        @test param.rtol == 1e-8
        @test param.atol == 1e-6
        @test param.zero_tol == 1e-5
        @test param.w0 == 1.0  # Default value
        
        # Test w0 parameter
        param_w0 = UnifiedOracleParam(w0 = 2.5)
        @test param_w0.w0 == 2.5
        
        # Test w0 validation - must be positive
        @test_throws ArgumentError UnifiedOracleParam(w0 = 0.0)
        @test_throws ArgumentError UnifiedOracleParam(w0 = -1.0)
        
        data = create_unified_test_data()
        
        function customize_sub_param!(model::Model, data::SimpleUnifiedTestData, scen_idx::Int; x)
            optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            I, J = data.n_facilities, data.n_customers
            @variable(model, y[1:I, 1:J] >= 0)
            
            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= x[i])
            
            return nothing
        end
        
        master = Master(data; customize = customize_master_unified!)
        oracle = UnifiedOracle(data, master; customize = customize_sub_param!, param = param)
        
        @test oracle.param.rtol == 1e-8
        @test oracle.param.zero_tol == 1e-5
        @test oracle.param.w0 == 1.0
    end
    
    @testset "UnifiedOracle with custom w0" begin
        # Test that UnifiedOracle works correctly with custom w0
        data = create_unified_test_data()
        mip_opt = solve_mip_unified_reference(data)
        
        function customize_sub_w0!(model::Model, data::SimpleUnifiedTestData, scen_idx::Int; x)
            optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            I, J = data.n_facilities, data.n_customers
            @variable(model, y[1:I, 1:J] >= 0)
            
            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= x[i])
            
            return nothing
        end
        
        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        
        # Test with w0 = 2.0
        param_w0 = UnifiedOracleParam(w0 = 2.0)
        master = Master(data; customize = customize_master_unified!)
        oracle = UnifiedOracle(data, master; customize = customize_sub_w0!, param = param_w0)
        env = BendersSeq(master, oracle; param = benders_param)
        solve!(env)
        
        @test env.termination_status == Optimal()
        @test isapprox(mip_opt, env.obj_value, atol=1e-4)
    end
    
    @testset "Error Handling" begin
        data = create_unified_test_data()

        # Note: QuadExpr constraint test is skipped because HiGHS doesn't support quadratic
        # constraints, so JuMP throws ErrorException before our code can process it.
        # This test would require a QP-capable solver like CPLEX or Gurobi.

        # Test 1: UnsupportedModelException for unsupported constraint set type (Interval)
        @testset "Interval constraint throws UnsupportedModelException" begin
            function customize_sub_interval_error!(model::Model, data::SimpleUnifiedTestData, scen_idx::Int; x)
                optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
                set_optimizer(model, optimizer)

                I, J = data.n_facilities, data.n_customers
                @variable(model, y[1:I, 1:J] >= 0)

                @objective(model, Min, sum(data.costs .* y))
                @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
                # Interval constraint - should trigger UnsupportedModelException
                @constraint(model, interval_con, 0 <= y[1,1] + x[1] <= 2)

                return nothing
            end

            master = Master(data; customize = customize_master_unified!)
            @test_throws UnsupportedModelException UnifiedOracle(data, master; customize = customize_sub_interval_error!)
        end

        # Test 2: ArgumentError for invalid w0
        @testset "Invalid w0 throws ArgumentError" begin
            @test_throws ArgumentError UnifiedOracleParam(w0 = 0.0)
            @test_throws ArgumentError UnifiedOracleParam(w0 = -1.0)
            @test_throws ArgumentError UnifiedOracleParam(w0 = -100.0)
        end
    end

    @testset "Reformulation Verification" begin
        # Test that the unified reformulation produces the correct model structure

        # Simple test problem with known structure
        struct ReformTestData <: AbstractData
            n::Int  # number of decision variables
        end

        function customize_master_reform!(model::Model, data::ReformTestData)
            optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
            set_optimizer(model, optimizer)

            @variable(model, x[1:data.n], Bin)
            @variable(model, t >= -1e6)
            @objective(model, Min, t)
            @constraint(model, sum(x) >= 1)

            return (x = x, ), t
        end

        @testset "σ variable structure" begin
            # Test that σ variable is added with correct bounds (σ >= 0)
            data = ReformTestData(2)

            function customize_sub_sigma!(model::Model, data::ReformTestData, scen_idx::Int; x)
                optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
                set_optimizer(model, optimizer)

                @variable(model, y >= 0)
                @objective(model, Min, y)
                @constraint(model, con1, y >= x[1] + x[2])

                return nothing
            end

            master = Master(data; customize = customize_master_reform!)
            oracle = UnifiedOracle(data, master; customize = customize_sub_sigma!)

            # Find σ variable by name
            σ = variable_by_name(oracle.model, "σ")
            @test σ !== nothing
            @test has_lower_bound(σ)
            @test lower_bound(σ) == 0.0
            @test !has_upper_bound(σ)
        end

        @testset "Objective becomes min σ" begin
            data = ReformTestData(2)

            function customize_sub_obj!(model::Model, data::ReformTestData, scen_idx::Int; x)
                optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
                set_optimizer(model, optimizer)

                @variable(model, y >= 0)
                @objective(model, Min, 3.0 * y + 5.0)  # Original objective
                @constraint(model, con1, y >= x[1] + x[2])

                return nothing
            end

            master = Master(data; customize = customize_master_reform!)
            oracle = UnifiedOracle(data, master; customize = customize_sub_obj!)

            # New objective should be just σ
            σ = variable_by_name(oracle.model, "σ")
            obj = objective_function(oracle.model)
            @test obj == σ
            @test objective_sense(oracle.model) == MOI.MIN_SENSE
        end

        @testset "Fixing constraints structure" begin
            # Verify lb (x + σ >= value) and ub (x - σ <= value) constraints
            data = ReformTestData(2)

            function customize_sub_fix!(model::Model, data::ReformTestData, scen_idx::Int; x)
                optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
                set_optimizer(model, optimizer)

                @variable(model, y >= 0)
                @objective(model, Min, y)
                @constraint(model, con1, y >= x[1] + x[2])

                return nothing
            end

            master = Master(data; customize = customize_master_reform!)
            oracle = UnifiedOracle(data, master; customize = customize_sub_fix!)

            # Should have one lb and one ub constraint per decision variable
            @test length(oracle.fixing_lb_constraints) == data.n
            @test length(oracle.fixing_ub_constraints) == data.n

            σ = variable_by_name(oracle.model, "σ")

            # Verify lb constraints: x + σ >= 0 (initial RHS)
            for con in oracle.fixing_lb_constraints
                con_obj = constraint_object(con)
                @test con_obj.set isa MOI.GreaterThan
                # σ coefficient should be +1
                @test normalized_coefficient(con, σ) == 1.0
            end

            # Verify ub constraints: x - σ <= 0 (initial RHS)
            for con in oracle.fixing_ub_constraints
                con_obj = constraint_object(con)
                @test con_obj.set isa MOI.LessThan
                # σ coefficient should be -1
                @test normalized_coefficient(con, σ) == -1.0
            end
        end

        @testset "Objective constraint structure" begin
            # Verify original objective becomes constraint: -obj + w0*σ >= -t
            data = ReformTestData(2)
            w0 = 2.5

            function customize_sub_objcon!(model::Model, data::ReformTestData, scen_idx::Int; x)
                optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
                set_optimizer(model, optimizer)

                @variable(model, y >= 0)
                @objective(model, Min, 3.0 * y)  # Original objective
                @constraint(model, con1, y >= x[1] + x[2])

                return nothing
            end

            master = Master(data; customize = customize_master_reform!)
            param = UnifiedOracleParam(w0 = w0)
            oracle = UnifiedOracle(data, master; customize = customize_sub_objcon!, param = param)

            σ = variable_by_name(oracle.model, "σ")

            # objective_constraint: -original_obj + w0*σ >= -t
            con_obj = constraint_object(oracle.objective_constraint)
            @test con_obj.set isa MOI.GreaterThan

            # σ coefficient should be w0
            @test normalized_coefficient(oracle.objective_constraint, σ) == w0
        end

        @testset "Problem constraints with σ relaxation" begin
            # Test that problem constraints are correctly modified with ±σ
            data = ReformTestData(2)

            function customize_sub_relax!(model::Model, data::ReformTestData, scen_idx::Int; x)
                optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
                set_optimizer(model, optimizer)

                @variable(model, y >= 0)
                @objective(model, Min, y)
                # >= constraint should get +σ
                @constraint(model, geq_con, y + x[1] >= 5.0)
                # <= constraint should get -σ
                @constraint(model, leq_con, y + x[2] <= 10.0)

                return nothing
            end

            master = Master(data; customize = customize_master_reform!)
            oracle = UnifiedOracle(data, master; customize = customize_sub_relax!)

            σ = variable_by_name(oracle.model, "σ")

            # Find >= constraint and verify σ coefficient is +1
            geq_con = constraint_by_name(oracle.model, "geq_con")
            @test geq_con !== nothing
            @test normalized_coefficient(geq_con, σ) == 1.0

            # Find <= constraint and verify σ coefficient is -1
            leq_con = constraint_by_name(oracle.model, "leq_con")
            @test leq_con !== nothing
            @test normalized_coefficient(leq_con, σ) == -1.0
        end

        @testset "Equality constraints split into lb/ub" begin
            # Test that == constraints are split into >= and <= with ±σ
            data = ReformTestData(2)

            function customize_sub_eq!(model::Model, data::ReformTestData, scen_idx::Int; x)
                optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
                set_optimizer(model, optimizer)

                @variable(model, y >= 0)
                @objective(model, Min, y)
                # == constraint should be split
                @constraint(model, eq_con, y + x[1] == 5.0)

                return nothing
            end

            master = Master(data; customize = customize_master_reform!)
            oracle = UnifiedOracle(data, master; customize = customize_sub_eq!)

            σ = variable_by_name(oracle.model, "σ")

            # Original eq_con should be deleted, replaced by eq_con_lb and eq_con_ub
            @test constraint_by_name(oracle.model, "eq_con") === nothing

            # lb version should exist with σ coefficient +1
            eq_con_lb = constraint_by_name(oracle.model, "eq_con_lb")
            @test eq_con_lb !== nothing
            con_obj_lb = constraint_object(eq_con_lb)
            @test con_obj_lb.set isa MOI.GreaterThan
            @test normalized_coefficient(eq_con_lb, σ) == 1.0

            # ub version should exist with σ coefficient -1
            eq_con_ub = constraint_by_name(oracle.model, "eq_con_ub")
            @test eq_con_ub !== nothing
            con_obj_ub = constraint_object(eq_con_ub)
            @test con_obj_ub.set isa MOI.LessThan
            @test normalized_coefficient(eq_con_ub, σ) == -1.0
        end

        @testset "Constraints without decision variables not relaxed" begin
            # Constraints not containing x should NOT have σ added
            data = ReformTestData(2)

            function customize_sub_nox!(model::Model, data::ReformTestData, scen_idx::Int; x)
                optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
                set_optimizer(model, optimizer)

                @variable(model, y >= 0)
                @variable(model, z >= 0)
                @objective(model, Min, y + z)
                # Constraint with x - should be relaxed
                @constraint(model, with_x, y + x[1] >= 1.0)
                # Constraint without x - should NOT be relaxed
                @constraint(model, without_x, y + z >= 2.0)

                return nothing
            end

            master = Master(data; customize = customize_master_reform!)
            oracle = UnifiedOracle(data, master; customize = customize_sub_nox!)

            σ = variable_by_name(oracle.model, "σ")

            # Constraint with x should have σ
            with_x_con = constraint_by_name(oracle.model, "with_x")
            @test normalized_coefficient(with_x_con, σ) == 1.0

            # Constraint without x should NOT have σ
            without_x_con = constraint_by_name(oracle.model, "without_x")
            @test normalized_coefficient(without_x_con, σ) == 0.0
        end
    end
end
