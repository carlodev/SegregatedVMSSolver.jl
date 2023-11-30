using SegregatedVMSSolver
using MPI
using Test
using PartitionedArrays

include(joinpath("..","CommonsTests", "InitializeParamsTests.jl"))
include(joinpath("..","CommonsTests","AddNewTagsTests.jl"))
include(joinpath("..","CommonsTests","StabParamsTests.jl"))
include(joinpath("..","CommonsTests","LinearUtilitiesTests.jl"))
include(joinpath("..","CommonsTests","StabilizedEquationsTests.jl"))
include(joinpath("..","CommonsTests","MatrixCreationTests.jl"))
include(joinpath("..","CommonsTests","RestartTests.jl"))

function test_common_mpi()
    with_mpi() do distribute
        comm = MPI.COMM_WORLD
        ranks =  MPI.Comm_rank(comm)
        #To avoid multiple printing of the same line in parallel
        if ranks != 0
            redirect_stderr(devnull)
            redirect_stdout(devnull)
        end
        InitializeParamsTests.main(distribute)
        AddNewTagsTests.main(distribute)
        StabParamsTests.main(distribute)
        LinearUtilitiesTests.main(distribute)
        StabilizedEquationsTests.main(distribute)
        MatrixCreationTests.main(distribute)
        RestartTests.main(distribute)

    end
end

function testmpi()
    @testset "Commons MPI" begin
        test_common_mpi()
        @test true
    end
end


testmpi()