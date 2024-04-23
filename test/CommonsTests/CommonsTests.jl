module CommonsTests
using SegregatedVMSSolver
using Test
using PartitionedArrays
using MPI



include(joinpath("CreateProblemTests", "CreateProblemTests.jl"))

include("MatrixCreationTests.jl")
include("EquationsTests.jl")


function test_solve(backend)
  @testset "Test Solve Problem $(backend)" begin
    backend() do distribute
      EquationsTests.main(distribute)
      MatrixCreationTests.main(distribute)
    end
  end
end

function tests_common(backend)
  if backend == with_mpi
    comm = MPI.COMM_WORLD
    ranks = MPI.Comm_rank(comm)
    #To avoid multiple printing of the same line in parallel
    if ranks != 0
      redirect_stderr(devnull)
      redirect_stdout(devnull)
    end
  end

  CreateProblemTests.test_create_problem(backend)
  test_solve(backend)
end






# pdir = joinpath(@__DIR__,"..","..",".")

# procs = 4
# function run_driver(procs,file)
#     mpiexec() do cmd

#          run(`$cmd -n $procs $(Base.julia_cmd()) --project=$pdir $file`)

#       @test true

#   end
# end

# run_driver(procs, joinpath("test","TestMPI.jl"))


end