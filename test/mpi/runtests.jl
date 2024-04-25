using SegregatedVMSSolver
using Parameters,PartitionedArrays
using Revise
using Test
using MPI

@testset "MPI" begin include( joinpath("MPITests.jl")) end



function clean_directory()
    rm_folders = ["Initial_Conditions", "Results", "Results_vtu"]
    # rm_paths = joinpath.("..",rm_folders)
    rm.(rm_folders;force=true, recursive=true)
end

#Clean outputs files/dir after testing
clean_directory()