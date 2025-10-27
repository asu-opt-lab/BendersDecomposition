using BendersDecomposition
using Test

@testset "BendersDecomposition.jl" begin
    include("1_test_sequential_typical/runtest.jl")
    include("2_test_sequential_in_out_typical/runtest.jl")
    include("3_test_sequential_disjunctive/runtest.jl")
    include("4_test_callback_typical/runtest.jl")
    include("5_test_callback_disjunctive/runtest.jl")
end
