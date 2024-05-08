using SegregatedVMSSolver
using Parameters,PartitionedArrays
using Revise
using Test
using MPI

@testset "MPI" begin include( joinpath("MPITests.jl")) end
