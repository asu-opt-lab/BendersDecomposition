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
end
