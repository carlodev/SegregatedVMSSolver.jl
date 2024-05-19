module SegregatedVMSSolverTests

using Test


@testset "Sequential" begin  include("sequential/runtests.jl") end

@testset "MPI" begin  include("mpi/runtests.jl") end

@testset "Utils" begin  include("UtilsTests/runtests.jl") end



end