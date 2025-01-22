using SegregatedVMSSolver
using Parameters,PartitionedArrays
using Revise
using Test


@testset "Sequential" begin include( joinpath("SequentialTests.jl")) end

@testset "Utils" begin  include(joinpath("..","UtilsTests","runtests.jl")) end