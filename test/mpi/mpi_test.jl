module MPItest
using PartitionedArrays
using MPI
using Test

# @testset "Commons MPI" begin
#   include(joinpath("..", "CommonsTests", "CommonsTests.jl"))
#   CommonsTests.tests_common(with_mpi)
# end

@testset "Cases Tests MPI" begin
    include(joinpath("..", "TestsCases", "TestsCases.jl"))
    TestsCases.tests_cases(with_mpi)
end

end