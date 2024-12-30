using Test
using JuMP
using Gurobi
using BendersDecomposition

@testset "DCGLP Tests" begin
    # Create test data
    data = CFLPData(
        5, 10,  # 5 facilities, 10 customers
        rand(5) .* 200,     # Random capacities
        rand(10) .* 50,     # Random demands  
        rand(5) .* 100,     # Random fixed costs
        rand(5, 10) .* 10   # Random costs
    )

    @testset "Base DCGLP Model" begin
        # Test model creation with different norm types
        @testset "StandardNorm" begin
            cut_strategy = ClassicalCut()
            disjunction_system = DisjunctiveCut

(
                cut_strategy,
                StandardNorm(),
                PureDisjunctionCut(),
                false,
                false,
                false,
                false
            )
            dcglp = create_dcglp(data, disjunction_system)

            # Test basic structure
            @test dcglp isa DCGLP
            @test dcglp.model isa Model
            @test haskey(dcglp.γ_constraints

, :γ₀)
            @test haskey(dcglp.γ_constraints

, :γₓ)
            @test haskey(dcglp.γ_constraints

, :γₜ)

            # Test variables exist
            model = dcglp.model
            @test :τ in keys(JuMP.object_dictionary(model))
            @test :k₀ in keys(JuMP.object_dictionary(model))
            @test :kₓ in keys(JuMP.object_dictionary(model))
            @test :v₀ in keys(JuMP.object_dictionary(model))
            @test :vₓ in keys(JuMP.object_dictionary(model))
            @test :kₜ in keys(JuMP.object_dictionary(model))
            @test :vₜ in keys(JuMP.object_dictionary(model))
        end

        @testset "LNorm Types" begin
            for norm_type in [L1Norm(), L2Norm(), LInfNorm()]
                @testset "$(typeof(norm_type))" begin
                    cut_strategy = ClassicalCut()
                    disjunction_system = DisjunctiveCut

(
                        cut_strategy,
                        norm_type,
                        PureDisjunctionCut(),
                        false,
                        false,
                        false,
                        false
                    )
                    dcglp = create_dcglp(data, disjunction_system)

                    # Test basic structure
                    @test dcglp isa DCGLP
                    @test dcglp.model isa Model

                    # Test additional LNorm variables
                    model = dcglp.model
                    @test :sx in keys(JuMP.object_dictionary(model))
                    @test :st in keys(JuMP.object_dictionary(model))

                    # Test cone constraints
                    if norm_type isa L1Norm || norm_type isa LInfNorm
                        @test any(S <: MOI.NormInfinityCone for (F,S) in list_of_constraint_types(model))
                    elseif norm_type isa L2Norm
                        @test any(S <: MOI.SecondOrderCone for (F,S) in list_of_constraint_types(model))
                    end
                end
            end
        end
    end
end
