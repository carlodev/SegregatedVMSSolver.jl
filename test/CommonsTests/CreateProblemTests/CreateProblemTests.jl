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


function main(distribute)
    @testset "Create Problem" begin
            AddNewTagsTests.main(distribute)
            RestartTests.main(distribute)
            ModelCreationTests.main(distribute)
            BoundaryConditionsTests.main(distribute)
            SpaceConditionsTests.main(distribute)
            InitialConditionsTests.main(distribute)     
    end
end




end ## end module