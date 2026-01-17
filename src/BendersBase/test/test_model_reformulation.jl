"""
UnifiedOracle2 Test Suite

Each test validates the following:
1. The objective function is the scalar variable z.
2. All reformulated constraints are expressed using MOI.GreaterThan{Float64}.
3. The variable z appears in the epigraph constraint and in every reformulated constraint that involves a master variable, and does not appear in any reformulated constraint that does not involve a master variable.


Tests various scenarios for different types of master variable container:
1. Vector{VariableRef} + master variable not named as x
2. Matrix{VariableRef}
3. DenseAxisArray{VariableRef}
4. Error handling: ArgumentError

"""

using Test
using BendersBase
using JuMP
using CPLEX

function valdiate_expr_type(model)
    for (t1, t2) in list_of_constraint_types(model)
        return (t1 == AffExpr && t2 != MOI.GreaterThan{Float64} ? false : true)
    end
end

function validate_slack_in_cons(model,x)
    z = model[:z]
    for con in all_constraints(model, AffExpr, MOI.GreaterThan{Float64})
        if con == model[:epigraph]
            lhs = JuMP.constraint_object(con).func
            z ∉ keys(lhs.terms) && (throw(UndefError("Slack z does not appear in epigraph constraint")); return false)
        else
            lhs= JuMP.constraint_object(con).func
            if !isempty(intersect(name.(keys(lhs.terms)), name.(x)))
                z ∉ keys(lhs.terms) && (throw(UndefError("Slack z does not appear in the constraint that master variable appears")); return false)
            else
                z ∈ keys(lhs.terms) && (throw(UndefError("Slack z appears in the constraint that master variable does not appear")); return false)
            end
        end
    end
    return true
end

# =============================================================================
# Test Scenarios
# =============================================================================

@testset verbose = true "UnifiedOracle Test Suite" begin
    @testset "Scenario 1: Vector{VariableRef}" begin
        # Tests master variable container is Vector{VariableRef}
        struct TestDataVar <: AbstractData
            n_facilities::Int64
            n_customers::Int64
            trans::Matrix{Float64}
            fixed::Vector{Float64}
        end

        I, J = 4, 5
        fixed = [8.0, 5.0, 9.0, 7.0]
        trans = [
                3.0 2.0 4.0 5.0 3.0;
                2.0 3.0 2.0 4.0 2.0;
                6.0 5.0 3.0 2.0 4.0;
                4.0 3.0 5.0 2.0 3.0
                ]
       data = TestDataVar(I, J, trans, fixed)

        function mip(data::AbstractData)
            I, J = data.n_facilities, data.n_customers

            model = Model()
            optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)

            @variable(model, x[1:I], Bin)
            @variable(model, y[1:I, 1:J] >= 0)
            
            @constraint(model, assign[j=1:J], sum(y[i,j] for i=1:I) == 1)
            @constraint(model, open_req[i=1:I, j=1:J], y[i,j] <= x[i])

            @objective(model, Min, sum(data.trans[i,j]*y[i,j] for i=1:I,j=1:J) + sum(data.fixed[i]*x[i] for i=1:I))
            optimize!(model)

            return objective_value(model), model
        end

        function customize_master_model!(model::Model, data::AbstractData)

            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)

            @variable(model, u[1:I], Bin)
            @variable(model, t >= -1e6)

            @objective(model, Min, sum(sum(data.fixed[i]*u[i] for i=1:I)) + t)
            
            return (u = u, ), t
        end
        
        function customize_sub_model!(model::Model, data::AbstractData, scen_idx::Int; u) 

            optimizer = optimizer_with_attributes(
                CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)
        
            @variable(model, y[1:I, 1:J] >= 0)

            @constraint(model, assign[j=1:J], sum(y[i,j] for i=1:I) == 1)
            @constraint(model, open_req[i=1:I, j=1:J], y[i,j] <= u[i])

            @objective(model, Min, sum(data.trans[i,j]*y[i,j] for i=1:I,j=1:J))
        end

        master = Master(data; customize = customize_master_model!)
        x_copy = copy_variables!(master.model, master.x_tuple); x = var_from_tuple(x_copy)
        oracle = UnifiedOracle(data, master; customize = customize_sub_model!); model = oracle.model

        @test objective_function(model) == model[:z] 
        @test valdiate_expr_type(model)
        @test validate_slack_in_cons(model, x)

        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        env = BendersSeq(master, oracle; param = benders_param)
        solve!(env)

        mip_obj, mip_model = mip(data)
        @test isapprox(mip_obj, objective_value(master.model); atol=1e-6)
    end
    @testset "Scenario 2: Matrix{VariableRef}" begin
        # Tests master variable container is Matrix{VariableRef}
        struct TestDataMat <: AbstractData
            n_facilities::Int64
            n_customers::Int64
            n_capacity_level::Int64
            trans::Matrix{Float64}
            fixed::Vector{Tuple{Float64,Float64}}
            capacity::Vector{Float64}
        end

        I, J, K = 3, 6, 2

        trans = [
            3.0 2.0 4.0 5.0 3.0 2.0;
            2.0 3.0 2.0 4.0 2.0 3.0;
            6.0 5.0 3.0 2.0 4.0 5.0
        ]
        fixed = [(6.0, 10.0), (5.0, 9.0), (8.0, 12.0)]
        capacity = [2.0, 4.0]
        data = TestDataMat(I, J, K, trans, fixed, capacity)

        function mip(data::AbstractData)
            I, J, K = data.n_facilities, data.n_customers, data.n_capacity_level

            model = Model()
            optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)

            @variable(model, x[1:I, 1:K], Bin)
            @variable(model, y[1:I, 1:J] >= 0)

            @constraint(model, assign[j=1:J], sum(y[i,j] for i=1:I) == 1)
            @constraint(model, one_level[i=1:I], sum(x[i,k] for k=1:K) <= 1)
            @constraint(model, capacity[i=1:I], sum(y[i,j] for j=1:J) <= sum(data.capacity[k]*x[i,k] for k=1:K))
            @constraint(model, open_req[i=1:I, j=1:J], y[i,j] <= sum(x[i,k] for k=1:K))

            @objective(model, Min,
                sum(data.trans[i,j]*y[i,j] for i=1:I,j=1:J) +
                sum(data.fixed[i][k] * x[i,k] for i=1:I, k=1:K)
            )

            optimize!(model)
            return objective_value(model), model
        end

        function customize_master_model!(model::Model, data::AbstractData)

            I, K = data.n_facilities, data.n_capacity_level

            optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)

            @variable(model, x[1:I, 1:K], Bin)
            @variable(model, t >= -1e6)
 
            @constraint(model, one_level[i=1:I], sum(x[i,k] for k=1:K) <= 1)

            @objective(model, Min, sum(data.fixed[i][k] * x[i,k] for i=1:I, k=1:K) + t)
            
            return (x = x, ), t
        end
        
        function customize_sub_model!(model::Model, data::AbstractData, scen_idx::Int; x) 
            I, J, K = data.n_facilities, data.n_customers, data.n_capacity_level

            optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)
        
            @variable(model, y[1:I, 1:J] >= 0)

            @constraint(model, assign[j=1:J], sum(y[i,j] for i=1:I) == 1)
            @constraint(model, capacity[i=1:I], sum(y[i,j] for j=1:J) <= sum(data.capacity[k]*x[i,k] for k=1:K))
            @constraint(model, open_req[i=1:I, j=1:J], y[i,j] <= sum(x[i,k] for k=1:K))

            @objective(model, Min, sum(data.trans[i,j]*y[i,j] for i=1:I,j=1:J))
        end

        master = Master(data; customize = customize_master_model!)
        x_copy = copy_variables!(master.model, master.x_tuple); x = var_from_tuple(x_copy)
        oracle = UnifiedOracle(data, master; customize = customize_sub_model!); model = oracle.model
        
        @test objective_function(model) == model[:z] 
        @test valdiate_expr_type(model)
        @test validate_slack_in_cons(model, x)

        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        env = BendersSeq(master, oracle; param = benders_param)
        solve!(env)

        mip_obj, mip_model = mip(data)
        @test isapprox(mip_obj, objective_value(master.model); atol=1e-6)
    end
    @testset "Scenario 3: DenseAxisArray{VariableRef}" begin
        # Tests master variable container is DenseAxisArray{VariableRef}
        struct TestDataDense <: AbstractData
            n_facilities::Int
            n_customers::Int
            modes::Vector{Symbol}
            trans::Matrix{Float64}
            fixed::Dict{Symbol, Vector{Float64}}
            capacity::Dict{Symbol, Float64}
        end

        I, J = 3, 5
        modes = [:small, :large]
        trans = [
            3.0 2.0 4.0 5.0 3.0;
            2.0 3.0 2.0 4.0 2.0;
            6.0 5.0 3.0 2.0 4.0
        ]
        fixed = Dict(:small => [6.0, 5.0, 8.0], :large => [10.0, 9.0, 12.0])
        capacity = Dict(:small => 2.0, :large => 4.0)

        data = TestDataDense(I, J, modes, trans, fixed, capacity)

        function mip(data::AbstractData)
            I, J = data.n_facilities, data.n_customers

            model = Model()
            optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)

            @variable(model, x[i=1:I, m in modes], Bin)
            @variable(model, y[1:I, 1:J] >= 0)

            @constraint(model, assign[j=1:J], sum(y[i,j] for i=1:I) == 1)
            @constraint(model, one_mode[i=1:I], sum(x[i,m] for m in data.modes) <= 1)
            @constraint(model, capacity[i=1:I], sum(y[i,j] for j=1:J) <= sum(data.capacity[m]*x[i,m] for m in modes))
            @constraint(model, open_req[i=1:I, j=1:J], y[i,j] <= sum(x[i,m] for m in data.modes))

            @objective(model, Min,
                sum(data.trans[i,j]*y[i,j] for i=1:I, j=1:J) +
                sum(data.fixed[m][i] * x[i,m] for i=1:I, m in data.modes)
            )
            optimize!(model)
            return objective_value(model), model
        end

        function customize_master_model!(model::Model, data::AbstractData)
            I = data.n_facilities

            optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)

            @variable(model, x[i=1:I, m in data.modes], Bin)
            @variable(model, t >= -1e6)

            @constraint(model, one_mode[i=1:I], sum(x[i,m] for m in data.modes) <= 1)
            
            @objective(model, Min, sum(data.fixed[m][i] * x[i,m] for i=1:I, m in data.modes) + t)
            
            return (x = x, ), t
        end
        
        function customize_sub_model!(model::Model, data::AbstractData, scen_idx::Int; x) 
            I, J = data.n_facilities, data.n_customers

            optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 1, MOI.Silent() => true)
            set_optimizer(model, optimizer)
        
            @variable(model, y[1:I, 1:J] >= 0)
            @constraint(model, assign[j=1:J], sum(y[i,j] for i=1:I) == 1)
            @constraint(model, capacity[i=1:I], sum(y[i,j] for j=1:J) <= sum(data.capacity[m]*x[i,m] for m in data.modes))
            @constraint(model, open_req[i=1:I, j=1:J], y[i,j] <= sum(x[i,m] for m in data.modes))

            @objective(model, Min, sum(data.trans[i,j]*y[i,j] for i=1:I, j=1:J))
        end

        master = Master(data; customize = customize_master_model!)
        x_copy = copy_variables!(master.model, master.x_tuple); x = var_from_tuple(x_copy)
        oracle = UnifiedOracle(data, master; customize = customize_sub_model!); model = oracle.model

        @test objective_function(model) == model[:z] 
        @test valdiate_expr_type(model)
        @test validate_slack_in_cons(model, x)

        benders_param = BendersSeqParam(time_limit = 60.0, gap_tolerance = 1e-6, verbose = false)
        env = BendersSeq(master, oracle; param = benders_param)
        solve!(env)

        mip_obj, mip_model = mip(data)
        @test isapprox(mip_obj, objective_value(master.model); atol=1e-6)
    end
    @testset "Scenario 4: Error Handling" begin
        # Tests when quadratic constraint is provided
        struct EmptyData <: AbstractData end
        data = EmptyData()

        function customize_master_model!(model::Model, data::EmptyData)
            @variable(model, u[1:2], Bin)
            @variable(model, t >= -1e6)
            return (u = u, ), t
        end
        
        function customize_sub_model1!(model::Model, data::EmptyData, scen_idx::Int; u) 
            @variable(model, y[1:2] >= 0)
            @constraint(model, y[1]^2 .<= u[1])
        end

        master = Master(data; customize = customize_master_model!)

        quad_model = Model()
        x_copy = copy_variables!(quad_model, master.x_tuple)
        customize_sub_model1!(quad_model, data, 1; x_copy...)
        x = var_from_tuple(x_copy)
        
        @test_throws ArgumentError model_reformulation!(quad_model, 1.0; x)
    end
end