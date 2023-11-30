using SegregatedVMSSolver
using Parameters,PartitionedArrays
using Revise
using Test
using MPI

@testset "Sequential" begin include( joinpath("sequential","SequentialTests.jl")) end
@testset "MPI" begin include( joinpath("mpi","MPITests.jl")) end

@testset "Utils" begin include(joinpath("UtilsTests","UtilsTests.jl")) end
