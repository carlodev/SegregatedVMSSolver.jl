# using SegregatedVMSSolver
# using Parameters,PartitionedArrays
# using Revise
# using Test
# using MPI

# @testset "Sequential" begin include( joinpath("sequential","SequentialTests.jl")) end
# @testset "MPI" begin include( joinpath("mpi","MPITests.jl")) end



# @testset "Utils" begin include(joinpath("UtilsTests","UtilsTests.jl")) end

function clean_directory()
    rm_folders = ["Initial_Conditions", "Results", "Results_vtu"]
    rm.(rm_folders;force=true, recursive=true)
end

#Clean outputs files/dir after testing
clean_directory()