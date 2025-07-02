using SegregatedVMSSolver
using Parameters,PartitionedArrays
using Revise
using Test


@testset "Sequential" begin include( joinpath("SequentialTests.jl")) end

