module CreateProblemTests

using SegregatedVMSSolver
using Test
using PartitionedArrays
using MPI


include("AddNewTagsTests.jl")
include("RestartTests.jl")

include("ModelCreationTests.jl")
include("BoundaryConditionsTests.jl")
include("SpaceConditionsTests.jl")
include("InitialConditionsTests.jl")


function test_create_problem(backend)
    @testset "Create Problem $(backend)" begin
        backend() do distribute
            AddNewTagsTests.main(distribute)
            RestartTests.main(distribute)
            ModelCreationTests.main(distribute)
            BoundaryConditionsTests.main(distribute)
            SpaceConditionsTests.main(distribute)
            InitialConditionsTests.main(distribute)
        end
    end
end


end ## end module