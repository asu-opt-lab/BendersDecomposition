"""
GBC (Generalized Bound Constraints) Test Suite

Tests various scenarios for GBC support:
1. Basic upper bound constraint: y <= x
2. Lower bound constraint: y >= 0.1*x
3. Fixed bound constraint: y == x
4. Scaled single variable: y <= 2*x
5. Multi-variable affine combination: y <= 0.5*x[1] + 0.5*x[2]
6. Partial y constraints: only subset of y has GBC
7. Mixed bound types: some upper, some lower, some fixed
8. No GBC (nothing return)
9. Mixed RHS types: mixing VariableRef and AffExpr for optimization testing
10. Error handling: DimensionMismatch, ArgumentError, Warning on empty LHS

Uses a simple transportation-like problem for testing.
"""

using BendersX
using Test
using JuMP
using CPLEX

# =============================================================================
# Test Problem: Simple Assignment Problem
# Minimize: sum(c[i,j] * y[i,j])
# Subject to: sum(y[:,j]) == 1 for all j (demand satisfied)
#            y[i,j] >= 0
#            y[i,j] <= x[i] for all i,j (facility open constraint - GBC)
# =============================================================================

struct SimpleAssignmentData <: AbstractData
    n_facilities::Int
    n_customers::Int
    costs::Matrix{Float64}
end

function create_test_data(n_facilities::Int=3, n_customers::Int=4)
    costs = rand(n_facilities, n_customers) .* 10
    return SimpleAssignmentData(n_facilities, n_customers, costs)
end

function customize_master_simple!(model::Model, data::SimpleAssignmentData)
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
    set_optimizer(model, optimizer)
    
    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t >= -1e6)
    @objective(model, Min, t)
    
    # At least one facility must be open
    @constraint(model, sum(x) >= 1)
    
    return (x = x, ), t
end

function solve_mip_reference(data::SimpleAssignmentData)
    model = Model()
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
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

@testset verbose = true "GBC Test Suite" begin
    
    @testset "Scenario 1: Basic Upper Bound" begin
        # Tests basic upper bound constraint y <= x
        data = create_test_data()
        mip_opt = solve_mip_reference(data)
        
        function customize_sub_upper!(model::Model, data::SimpleAssignmentData, scen_idx::Int; x)
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            I, J = data.n_facilities, data.n_customers
            @variable(model, y[1:I, 1:J] >= 0)
            
            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            
            # Unified format: gbc_rhs is Vector{AffExpr}
            gbc_lhs = vec(y)
            gbc_rhs = [AffExpr(0.0, x[i] => 1.0) for i in 1:I for j in 1:J]
            gbc_sense = fill(UpperBound, I*J)
            
            return gbc_lhs, gbc_rhs, gbc_sense
        end
        
        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        master = Master(data; customize = customize_master_simple!)
        oracle = ClassicalOracle(data, master; customize = customize_sub_upper!)
        env = BendersSeq(master, oracle; param = benders_param)
        solve!(env)
        
        @test env.termination_status == Optimal()
        @test isapprox(mip_opt, env.obj_value, atol=1e-4)
    end
    
    @testset "Scenario 2: Lower Bound Constraint" begin
        # Tests y >= x constraint (facility must serve if open)
        # Modified problem: if facility is open, it must serve at least some amount
        data = create_test_data()

        function solve_mip_scenario_3(data::SimpleAssignmentData)
            model = Model()
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)

            I, J = data.n_facilities, data.n_customers
            @variable(model, x[1:I], Bin)
            @variable(model, y[1:I, 1:J] >= 0)
            @variable(model, z[1:I] >= 0)  # Total service from facility i

            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= 1)
            @constraint(model, total_service[i in 1:I], z[i] == sum(y[i,:]))
            @constraint(model, min_service[i in 1:I], z[i] >= 0.1 * x[i])  # GBC: Lower bound
            @constraint(model, sum(x) >= 1)

            optimize!(model)
            return objective_value(model)
        end

        mip_opt = solve_mip_scenario_3(data)
        
        function customize_sub_lower!(model::Model, data::SimpleAssignmentData, scen_idx::Int; x)
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            I, J = data.n_facilities, data.n_customers
            @variable(model, y[1:I, 1:J] >= 0)
            @variable(model, z[1:I] >= 0)  # Total service from facility i
            
            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= 1)  # Upper bound
            @constraint(model, total_service[i in 1:I], z[i] == sum(y[i,:]))
            
            # Lower bound GBC: z[i] >= 0.1 * x[i] (if open, serve at least 10%)
            gbc_lhs = collect(z)
            gbc_rhs = [0.1 * x[i] for i in 1:I]
            gbc_sense = fill(LowerBound, I)
            
            return gbc_lhs, gbc_rhs, gbc_sense
        end
        
        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        master = Master(data; customize = customize_master_simple!)
        oracle = ClassicalOracle(data, master; customize = customize_sub_lower!)
        env = BendersSeq(master, oracle; param = benders_param)
        solve!(env)

        @test env.termination_status == Optimal()
        @test isapprox(mip_opt, env.obj_value, atol=1e-4)
    end

    @testset "Scenario 3: Fixed Bound Constraint" begin
        # Tests y == x constraint
        data = create_test_data()

        function solve_mip_scenario_4(data::SimpleAssignmentData)
            model = Model()
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)

            I, J = data.n_facilities, data.n_customers
            @variable(model, x[1:I], Bin)
            @variable(model, y[1:I, 1:J] >= 0)
            @variable(model, w[1:I] >= 0)  # Auxiliary: exactly equals x

            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= w[i])
            @constraint(model, fixed_capacity[i in 1:I], w[i] == x[i])  # GBC: Fixed bound
            @constraint(model, sum(x) >= 1)

            optimize!(model)
            return objective_value(model)
        end

        mip_opt = solve_mip_scenario_4(data)
        
        function customize_sub_fixed!(model::Model, data::SimpleAssignmentData, scen_idx::Int; x)
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            I, J = data.n_facilities, data.n_customers
            @variable(model, y[1:I, 1:J] >= 0)
            @variable(model, w[1:I] >= 0)  # Auxiliary: exactly equals x
            
            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= w[i])
            
            # Fixed bound GBC: w[i] == x[i]
            gbc_lhs = collect(w)
            gbc_rhs = [1.0 * x[i] for i in 1:I]
            gbc_sense = fill(FixedBound, I)
            
            return gbc_lhs, gbc_rhs, gbc_sense
        end
        
        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        master = Master(data; customize = customize_master_simple!)
        oracle = ClassicalOracle(data, master; customize = customize_sub_fixed!)
        env = BendersSeq(master, oracle; param = benders_param)
        solve!(env)

        @test env.termination_status == Optimal()
        @test isapprox(mip_opt, env.obj_value, atol=1e-4)
    end

    @testset "Scenario 4: Scaled Single Variable" begin
        # Tests y <= 2*x constraint
        data = create_test_data()

        function solve_mip_scenario_5(data::SimpleAssignmentData)
            model = Model()
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)

            I, J = data.n_facilities, data.n_customers
            @variable(model, x[1:I], Bin)
            @variable(model, y[1:I, 1:J] >= 0)

            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= 2.0 * x[i])  # GBC: Scaled (2.0)
            @constraint(model, sum(x) >= 1)

            optimize!(model)
            return objective_value(model)
        end

        mip_opt = solve_mip_scenario_5(data)
        
        function customize_sub_scaled!(model::Model, data::SimpleAssignmentData, scen_idx::Int; x)
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            I, J = data.n_facilities, data.n_customers
            @variable(model, y[1:I, 1:J] >= 0)
            
            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            
            # Scaled GBC: y[i,j] <= 2*x[i] (allow up to 2 units per facility)
            gbc_lhs = vec(y)
            gbc_rhs = [2.0 * x[i] for i in 1:I for j in 1:J]
            gbc_sense = fill(UpperBound, I*J)
            
            return gbc_lhs, gbc_rhs, gbc_sense
        end
        
        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        master = Master(data; customize = customize_master_simple!)
        oracle = ClassicalOracle(data, master; customize = customize_sub_scaled!)
        env = BendersSeq(master, oracle; param = benders_param)
        solve!(env)

        @test env.termination_status == Optimal()
        @test isapprox(mip_opt, env.obj_value, atol=1e-4)
    end

    @testset "Scenario 5: Multi-Variable Affine Combination" begin
        # Tests y <= x[1] + x[2] constraint (depends on multiple x)
        data = create_test_data(4, 3)  # 4 facilities, 3 customers

        function solve_mip_scenario_6(data::SimpleAssignmentData)
            model = Model()
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)

            I, J = data.n_facilities, data.n_customers
            @variable(model, x[1:I], Bin)
            @variable(model, y[1:J] >= 0)  # Total service to each customer

            @objective(model, Min, sum(data.costs[1,:] .* y))  # Use costs from facility 1
            @constraint(model, demand[j in 1:J], y[j] <= 1)
            @constraint(model, total, sum(y) >= 1)
            @constraint(model, affine_gbc[j in 1:J], y[j] <= 0.5*x[1] + 0.5*x[2])  # GBC: Affine combination
            @constraint(model, sum(x) >= 2)  # At least 2 facilities

            optimize!(model)
            return objective_value(model)
        end

        mip_opt = solve_mip_scenario_6(data)

        function customize_master_pairs!(model::Model, data::SimpleAssignmentData)
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            @variable(model, x[1:data.n_facilities], Bin)
            @variable(model, t >= -1e6)
            @objective(model, Min, t)
            @constraint(model, sum(x) >= 2)  # At least 2 facilities
            
            return (x = x, ), t
        end
        
        function customize_sub_affine!(model::Model, data::SimpleAssignmentData, scen_idx::Int; x)
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            I, J = data.n_facilities, data.n_customers
            @variable(model, y[1:J] >= 0)  # Total service to each customer
            
            @objective(model, Min, sum(data.costs[1,:] .* y))  # Use costs from facility 1
            @constraint(model, demand[j in 1:J], y[j] <= 1)
            @constraint(model, total, sum(y) >= 1)
            
            # Affine GBC: y[j] <= x[1] + x[2] for each customer j
            # This means customer j can be served only if facilities 1 OR 2 are open
            gbc_lhs = collect(y)
            gbc_rhs = [0.5*x[1] + 0.5*x[2] for _ in 1:J]  # y[j] <= 0.5*x[1] + 0.5*x[2]
            gbc_sense = fill(UpperBound, J)
            
            return gbc_lhs, gbc_rhs, gbc_sense
        end
        
        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        master = Master(data; customize = customize_master_pairs!)
        oracle = ClassicalOracle(data, master; customize = customize_sub_affine!)
        env = BendersSeq(master, oracle; param = benders_param)
        solve!(env)

        @test env.termination_status == Optimal()
        @test isapprox(mip_opt, env.obj_value, atol=1e-4)
    end

    @testset "Scenario 6: Partial Y Constraints" begin
        # Tests case where only subset of y has GBC constraints
        data = create_test_data()

        function solve_mip_scenario_7(data::SimpleAssignmentData)
            model = Model()
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)

            I, J = data.n_facilities, data.n_customers
            @variable(model, x[1:I], Bin)
            @variable(model, y[1:I, 1:J] >= 0)

            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            # GBC: Only first facility has constraint
            @constraint(model, facility_1_open[j in 1:J], y[1,j] <= x[1])
            # No constraints on y[2:I, :] with respect to x
            @constraint(model, sum(x) >= 1)

            optimize!(model)
            return objective_value(model)
        end

        mip_opt = solve_mip_scenario_7(data)
        
        function customize_sub_partial!(model::Model, data::SimpleAssignmentData, scen_idx::Int; x)
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            I, J = data.n_facilities, data.n_customers
            @variable(model, y[1:I, 1:J] >= 0)
            
            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            
            # Only first facility has GBC constraint
            # y[1,j] <= x[1] for all j
            gbc_lhs = [y[1,j] for j in 1:J]
            gbc_rhs = [1.0 * x[1] for _ in 1:J]
            gbc_sense = fill(UpperBound, J)
            
            return gbc_lhs, gbc_rhs, gbc_sense
        end
        
        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        master = Master(data; customize = customize_master_simple!)
        oracle = ClassicalOracle(data, master; customize = customize_sub_partial!)
        env = BendersSeq(master, oracle; param = benders_param)
        solve!(env)

        @test env.termination_status == Optimal()
        @test isapprox(mip_opt, env.obj_value, atol=1e-4)
    end

    @testset "Scenario 7: Mixed Bound Types" begin
        # Tests case with different bound types for different variables
        data = create_test_data()

        function solve_mip_scenario_8(data::SimpleAssignmentData)
            model = Model()
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)

            I, J = data.n_facilities, data.n_customers
            @variable(model, x[1:I], Bin)
            @variable(model, y[1:I, 1:J] >= 0)
            @variable(model, z[1:I] >= 0)  # Total service
            @variable(model, w[1:I] >= 0)  # Fixed capacity

            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= w[i])
            @constraint(model, total_service[i in 1:I], z[i] == sum(y[i,:]))
            # Mixed GBC:
            @constraint(model, fixed_capacity[i in 1:I], w[i] == x[i])  # Fixed bound
            @constraint(model, min_service[i in 1:I], z[i] >= 0.01 * x[i])  # Lower bound
            @constraint(model, sum(x) >= 1)

            optimize!(model)
            return objective_value(model)
        end

        mip_opt = solve_mip_scenario_8(data)
        
        function customize_sub_mixed!(model::Model, data::SimpleAssignmentData, scen_idx::Int; x)
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            I, J = data.n_facilities, data.n_customers
            @variable(model, y[1:I, 1:J] >= 0)
            @variable(model, z[1:I] >= 0)  # Total service
            @variable(model, w[1:I] >= 0)  # Fixed capacity
            
            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= w[i])
            @constraint(model, total_service[i in 1:I], z[i] == sum(y[i,:]))
            
            # Mixed GBC:
            # - w[i] == x[i] (Fixed)
            # - z[i] >= 0.01*x[i] (Lower, if open serve at least 1%)
            gbc_lhs = vcat(collect(w), collect(z))
            gbc_rhs = vcat([1.0 * x[i] for i in 1:I], [0.01 * x[i] for i in 1:I])
            gbc_sense = vcat(fill(FixedBound, I), fill(LowerBound, I))
            
            return gbc_lhs, gbc_rhs, gbc_sense
        end
        
        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        master = Master(data; customize = customize_master_simple!)
        oracle = ClassicalOracle(data, master; customize = customize_sub_mixed!)
        env = BendersSeq(master, oracle; param = benders_param)
        solve!(env)

        @test env.termination_status == Optimal()
        @test isapprox(mip_opt, env.obj_value, atol=1e-4)
    end

    @testset "Scenario 8: No GBC (Nothing Return)" begin
        # Tests that code works when customize returns nothing
        data = create_test_data()
        mip_opt = solve_mip_reference(data)
        
        function customize_sub_no_gbc!(model::Model, data::SimpleAssignmentData, scen_idx::Int; x)
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            I, J = data.n_facilities, data.n_customers
            @variable(model, y[1:I, 1:J] >= 0)
            
            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            @constraint(model, facility_open[i in 1:I, j in 1:J], y[i,j] <= x[i])  # Regular constraint
            
            # No GBC return
            return nothing
        end
        
        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        master = Master(data; customize = customize_master_simple!)
        oracle = ClassicalOracle(data, master; customize = customize_sub_no_gbc!)
        env = BendersSeq(master, oracle; param = benders_param)
        solve!(env)

        @test env.termination_status == Optimal()
        @test isapprox(mip_opt, env.obj_value, atol=1e-4)
    end

    @testset "Scenario 9: Mixed RHS Types (VariableRef + AffExpr)" begin
        # Tests mixing VariableRef and AffExpr in gbc_rhs
        data = create_test_data()
        mip_opt = solve_mip_reference(data)
        
        function customize_sub_union!(model::Model, data::SimpleAssignmentData, scen_idx::Int; x)
            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)
            
            I, J = data.n_facilities, data.n_customers
            @variable(model, y[1:I, 1:J] >= 0)
            
            @objective(model, Min, sum(data.costs .* y))
            @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
            
            # Use mixed Vector{Union{VariableRef, AffExpr}}
            # First half: use VariableRef
            # Second half: use AffExpr (1.0 * x)
            gbc_lhs = vec(y)
            gbc_rhs = Vector{Union{VariableRef, AffExpr}}(undef, length(y))
            
            idx = 1
            mid = length(y) ÷ 2
            for i in 1:I, j in 1:J
                if idx <= mid
                     gbc_rhs[idx] = x[i]  # VariableRef
                else
                     gbc_rhs[idx] = 1.0 * x[i] + 0.0  # AffExpr
                end
                idx += 1
            end
            
            gbc_sense = fill(UpperBound, I*J)
            return gbc_lhs, gbc_rhs, gbc_sense
        end
        
        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        master = Master(data; customize = customize_master_simple!)
        oracle = ClassicalOracle(data, master; customize = customize_sub_union!)
        env = BendersSeq(master, oracle; param = benders_param)
        solve!(env)
        
        @test env.termination_status == Optimal()
        @test isapprox(mip_opt, env.obj_value, atol=1e-4)
    end
    
    @testset "Scenario 10: Error Handling" begin
        data = create_test_data()
        
        # 1. Test DimensionMismatch
        function customize_sub_dim_error!(model::Model, data::SimpleAssignmentData, scen_idx::Int; x)
            gbc_lhs = [VariableRef(model, MOI.VariableIndex(1))]
            gbc_rhs = union_type([x[1], x[1]]) # Length 2
            gbc_sense = [UpperBound]
            return gbc_lhs, gbc_rhs, gbc_sense
        end
        
        master = Master(data; customize = customize_master_simple!)
        # Helper function to supply Union type
        union_type(v) = Vector{Union{VariableRef, AffExpr}}(v)

        @test_throws DimensionMismatch ClassicalOracle(data, master; customize = customize_sub_dim_error!)
        
        # 2. Test ArgumentError (Invalid tuple length)
        function customize_sub_arg_error!(model::Model, data::SimpleAssignmentData, scen_idx::Int; x)
            return ([x[1]], [x[1]]) # Length 2, expected 3
        end
        @test_throws ArgumentError ClassicalOracle(data, master; customize = customize_sub_arg_error!)

         # 3. Test Empty LHS Warning
        function customize_sub_empty_warn!(model::Model, data::SimpleAssignmentData, scen_idx::Int; x)
            return VariableRef[], Union{VariableRef, AffExpr}[], GBCBoundType[]
        end
        
        # Capture warning
        @test_logs (:warn, "GBC tuple returned but gbc_lhs is empty. No GBC constraints will be applied.") begin
            ClassicalOracle(data, master; customize = customize_sub_empty_warn!)
        end
    end

end

