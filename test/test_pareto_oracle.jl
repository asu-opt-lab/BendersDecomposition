"""
ParetoOracle Test Suite

Tests the ParetoOracle implementation for Benders decomposition.
Compares results between ClassicalOracle and ParetoOracle to ensure consistency.

Uses a simple assignment problem structure for testing.
"""

using BendersX
using Test
using JuMP
using HiGHS

# =============================================================================
# Test Problem: Simple Assignment Problem
# =============================================================================

struct SimpleParetoTestData <: AbstractData
    n_facilities::Int
    n_customers::Int
    costs::Matrix{Float64}
end

function create_pareto_test_data(n_facilities::Int=3, n_customers::Int=4)
    # Use fixed seed for reproducibility
    costs = [1.0 2.0 3.0 4.0;
             5.0 1.0 2.0 3.0;
             4.0 3.0 1.0 2.0]
    return SimpleParetoTestData(n_facilities, n_customers, costs)
end

function customize_master_pareto!(model::Model, data::SimpleParetoTestData)
    optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
    set_optimizer(model, optimizer)
    
    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t >= -1e6)
    @objective(model, Min, t)
    
    # At least one facility must be open
    @constraint(model, sum(x) >= 1)
    
    return (x = x, ), t
end

function solve_mip_pareto_reference(data::SimpleParetoTestData)
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

function customize_sub_pareto!(model::Model, data::SimpleParetoTestData, scen_idx::Int; x)
    optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
    set_optimizer(model, optimizer)
    
    I, J = data.n_facilities, data.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    
    @objective(model, Min, sum(data.costs .* y))
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= x[i])
    
    return nothing
end

# =============================================================================
# Test Scenarios
# =============================================================================

@testset verbose = true "ParetoOracle Test Suite" begin
    
    @testset "ParetoOracleParam Construction" begin
        # Test that core_point is required
        @test_throws MethodError ParetoOracleParam()
        
        # Test that empty core_point throws error
        @test_throws ArgumentError ParetoOracleParam(Float64[])
        
        # Test valid construction
        param = ParetoOracleParam([0.5, 0.5, 0.5])
        @test param.core_point == [0.5, 0.5, 0.5]
        @test param.rtol == 1e-9
        @test param.atol == 0.0
        @test param.zero_tol == 1e-9
        @test param.λ == 1.0  # Default λ
        
        # Test with custom tolerances
        param2 = ParetoOracleParam([0.3, 0.7]; rtol = 1e-8, atol = 1e-6, zero_tol = 1e-5)
        @test param2.core_point == [0.3, 0.7]
        @test param2.rtol == 1e-8
        @test param2.atol == 1e-6
        @test param2.zero_tol == 1e-5
        @test param2.λ == 1.0  # Default λ
        
        # Test with custom λ
        param3 = ParetoOracleParam([0.5, 0.5]; λ = 0.8)
        @test param3.λ == 0.8
        
        # Test λ must be in [0, 1]
        @test_throws ArgumentError ParetoOracleParam([0.5]; λ = -0.1)
        @test_throws ArgumentError ParetoOracleParam([0.5]; λ = 1.5)
    end
    
    @testset "Basic ParetoOracle Construction" begin
        data = create_pareto_test_data()
        master = Master(data; customize = customize_master_pareto!)
        
        # core_point dimension must match decision variables
        param = ParetoOracleParam(fill(0.5, data.n_facilities))
        oracle = ParetoOracle(data, master, param; customize = customize_sub_pareto!)
        
        # Verify oracle has all required fields
        @test !isempty(oracle.fixed_x_constraints)
        @test !isempty(oracle.pareto_fixing_constraints)
        @test oracle.pareto_variable isa VariableRef
        @test oracle.model isa Model
        @test oracle.pareto_model isa Model
        @test oracle.param.core_point == fill(0.5, data.n_facilities)
    end
    
    @testset "ParetoOracle Dimension Mismatch" begin
        data = create_pareto_test_data()
        master = Master(data; customize = customize_master_pareto!)
        
        # Wrong dimension for core_point should throw
        param_wrong = ParetoOracleParam([0.5, 0.5])  # Only 2 elements, but n_facilities = 3
        @test_throws DimensionMismatch ParetoOracle(data, master, param_wrong; customize = customize_sub_pareto!)
    end
    
    @testset "ParetoOracle Unsupported Constraint Type" begin
        # Test that interval constraints throw UnsupportedModelException
        data = create_pareto_test_data()
        master = Master(data; customize = customize_master_pareto!)
        param = ParetoOracleParam(fill(0.5, data.n_facilities))
        
        # Custom subproblem with interval constraint (unsupported)
        function customize_sub_with_interval!(model::Model, data::SimpleParetoTestData, scen_idx::Int; x)
            optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            I, J = data.n_facilities, data.n_customers
            @variable(model, y[1:I, 1:J] >= 0)
            
            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= x[i])
            
            # Add interval constraint (not supported by ParetoOracle)
            @constraint(model, 0.5 <= sum(y) <= 10.0)
            
            return nothing
        end
        
        @test_throws UnsupportedModelException ParetoOracle(data, master, param; customize = customize_sub_with_interval!)
    end
    
    @testset "ParetoOracle vs ClassicalOracle - Basic Problem" begin
        # Compare results between ClassicalOracle and ParetoOracle
        data = create_pareto_test_data()
        mip_opt = solve_mip_pareto_reference(data)
        
        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        
        # Test ClassicalOracle
        master_classical = Master(data; customize = customize_master_pareto!)
        oracle_classical = ClassicalOracle(data, master_classical; customize = customize_sub_pareto!)
        env_classical = BendersSeq(master_classical, oracle_classical; param = benders_param)
        solve!(env_classical)
        
        @test env_classical.termination_status == Optimal()
        @test isapprox(mip_opt, env_classical.obj_value, atol=1e-4)
        
        # Test ParetoOracle
        master_pareto = Master(data; customize = customize_master_pareto!)
        pareto_param = ParetoOracleParam(fill(0.5, data.n_facilities))
        oracle_pareto = ParetoOracle(data, master_pareto, pareto_param; customize = customize_sub_pareto!)
        env_pareto = BendersSeq(master_pareto, oracle_pareto; param = benders_param)
        solve!(env_pareto)
        
        @test env_pareto.termination_status == Optimal()
        @test isapprox(mip_opt, env_pareto.obj_value, atol=1e-4)
        
        # Both should find the same solution
        @test isapprox(env_classical.obj_value, env_pareto.obj_value, atol=1e-4)
    end
    
    @testset "ParetoOracle with Different Core Points" begin
        # Note: Core point must be in relative interior of feasible region for x
        # For this problem, x ∈ [0,1]^3 with sum(x) >= 1
        # So core_point = [0.5, 0.5, 0.5] is a valid interior point (sum = 1.5 > 1)
        
        data = create_pareto_test_data()
        mip_opt = solve_mip_pareto_reference(data)
        
        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        
        # Test with core_point = [0.5, 0.5, 0.5] (interior point)
        master1 = Master(data; customize = customize_master_pareto!)
        param1 = ParetoOracleParam(fill(0.5, data.n_facilities))
        oracle1 = ParetoOracle(data, master1, param1; customize = customize_sub_pareto!)
        env1 = BendersSeq(master1, oracle1; param = benders_param)
        solve!(env1)
        
        @test env1.termination_status == Optimal()
        @test isapprox(mip_opt, env1.obj_value, atol=1e-4)
        
        # Test with core_point = [0.8, 0.8, 0.8] (another interior point)
        master2 = Master(data; customize = customize_master_pareto!)
        param2 = ParetoOracleParam(fill(0.8, data.n_facilities))
        oracle2 = ParetoOracle(data, master2, param2; customize = customize_sub_pareto!)
        env2 = BendersSeq(master2, oracle2; param = benders_param)
        solve!(env2)
        
        @test env2.termination_status == Optimal()
        @test isapprox(mip_opt, env2.obj_value, atol=1e-4)
        
        # Test with asymmetric interior core_point (sum = 1.5 > 1)
        master3 = Master(data; customize = customize_master_pareto!)
        param3 = ParetoOracleParam([0.4, 0.5, 0.6])
        oracle3 = ParetoOracle(data, master3, param3; customize = customize_sub_pareto!)
        env3 = BendersSeq(master3, oracle3; param = benders_param)
        solve!(env3)
        
        @test env3.termination_status == Optimal()
        @test isapprox(mip_opt, env3.obj_value, atol=1e-4)
    end
    
    @testset "ParetoOracle generate_cuts Basic" begin
        data = create_pareto_test_data()
        master = Master(data; customize = customize_master_pareto!)
        param = ParetoOracleParam(fill(0.5, data.n_facilities))
        oracle = ParetoOracle(data, master, param; customize = customize_sub_pareto!)
        
        # Test generate_cuts with a feasible point
        x_value = [1.0, 0.0, 0.0]  # Only facility 1 open
        t_value = [0.0]  # Initial t value
        
        is_in_L, hyperplanes, sub_obj_vals = generate_cuts(oracle, x_value, t_value)
        
        # Should return valid hyperplane
        @test length(hyperplanes) == 1
        @test length(sub_obj_vals) == 1
        @test sub_obj_vals[1] >= 0  # Objective should be non-negative
        
        # Hyperplane should have correct dimensions
        @test length(hyperplanes[1].a_x) == data.n_facilities
        @test length(hyperplanes[1].a_t) == 1
    end
    
    @testset "ParetoOracle Dynamic Core Point Update (λ parameter)" begin
        # Test that λ < 1 updates core_point after each generate_cuts call
        data = create_pareto_test_data()
        master = Master(data; customize = customize_master_pareto!)
        
        # Test with λ = 0.5 (core_point should move towards x_value)
        initial_core = fill(0.5, data.n_facilities)
        param = ParetoOracleParam(copy(initial_core); λ = 0.5)
        oracle = ParetoOracle(data, master, param; customize = customize_sub_pareto!)
        
        # Verify initial core_point
        @test oracle.param.core_point == initial_core
        
        # First call to generate_cuts
        x_value_1 = [1.0, 0.0, 0.0]
        t_value = [0.0]
        generate_cuts(oracle, x_value_1, t_value)
        
        # core_point should now be: 0.5 * [0.5,0.5,0.5] + 0.5 * [1.0,0.0,0.0] = [0.75, 0.25, 0.25]
        expected_core_1 = 0.5 .* initial_core .+ 0.5 .* x_value_1
        @test isapprox(oracle.param.core_point, expected_core_1, atol=1e-10)
        
        # Second call with different x_value
        x_value_2 = [0.0, 1.0, 0.0]
        generate_cuts(oracle, x_value_2, t_value)
        
        # core_point should now be: 0.5 * [0.75,0.25,0.25] + 0.5 * [0.0,1.0,0.0] = [0.375, 0.625, 0.125]
        expected_core_2 = 0.5 .* expected_core_1 .+ 0.5 .* x_value_2
        @test isapprox(oracle.param.core_point, expected_core_2, atol=1e-10)
    end
    
    @testset "ParetoOracle with λ = 1.0 (no update)" begin
        # Test that λ = 1.0 keeps core_point unchanged (classical behavior)
        data = create_pareto_test_data()
        master = Master(data; customize = customize_master_pareto!)
        
        initial_core = fill(0.5, data.n_facilities)
        param = ParetoOracleParam(copy(initial_core); λ = 1.0)  # Default: no update
        oracle = ParetoOracle(data, master, param; customize = customize_sub_pareto!)
        
        # Verify initial core_point
        @test oracle.param.core_point == initial_core
        
        # Call generate_cuts multiple times
        x_value = [1.0, 0.0, 0.0]
        t_value = [0.0]
        generate_cuts(oracle, x_value, t_value)
        generate_cuts(oracle, x_value, t_value)
        
        # core_point should remain unchanged with λ = 1.0
        @test oracle.param.core_point == initial_core
    end
    
    @testset "ParetoOracle with λ = 0.0 (full update)" begin
        # Test that λ = 0.0 replaces core_point with x_value entirely
        data = create_pareto_test_data()
        master = Master(data; customize = customize_master_pareto!)
        
        initial_core = fill(0.5, data.n_facilities)
        param = ParetoOracleParam(copy(initial_core); λ = 0.0)  # Full update
        oracle = ParetoOracle(data, master, param; customize = customize_sub_pareto!)
        
        # Call generate_cuts
        x_value = [1.0, 0.0, 0.0]
        t_value = [0.0]
        generate_cuts(oracle, x_value, t_value)
        
        # core_point should be exactly x_value
        @test oracle.param.core_point == x_value
    end
    
    @testset "ParetoOracle with Dynamic Core Point - Solution Quality" begin
        # Verify that dynamic core point still produces correct solutions
        data = create_pareto_test_data()
        mip_opt = solve_mip_pareto_reference(data)
        
        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        
        # Test with λ = 0.8 (gradual update)
        master = Master(data; customize = customize_master_pareto!)
        param = ParetoOracleParam(fill(0.5, data.n_facilities); λ = 0.8)
        oracle = ParetoOracle(data, master, param; customize = customize_sub_pareto!)
        env = BendersSeq(master, oracle; param = benders_param)
        solve!(env)
        
        @test env.termination_status == Optimal()
        @test isapprox(mip_opt, env.obj_value, atol=1e-4)
    end
end
