using Test
using PartitionedArrays

@testset "Cases Tests Sequential" begin
    include(joinpath("..", "TestsCases", "TestsCases.jl"))
    TestsCases.tests_cases(with_debug)
end
