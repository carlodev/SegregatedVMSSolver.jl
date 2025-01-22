module SequentialTests
using Test
using PartitionedArrays

@testset "Commons Sequential" begin
    include(joinpath("..", "CommonsTests", "CommonsTests.jl"))
    CommonsTests.tests_common(with_debug)
end

@testset "Cases Tests Sequential" begin
    include(joinpath("..", "TestsCases", "TestsCases.jl"))
    TestsCases.tests_cases(with_debug)
end


end

