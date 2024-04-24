module TestsCases

using PartitionedArrays
using SegregatedVMSSolver
using Test
using MPI
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions

include("AirfoilTest.jl")
include("CylinderTest.jl")
include("LidDrivenTest.jl")
include("TaylorGreenTest.jl")

function tests_cases(backend)
  backend() do distribute
    if backend == with_mpi
        comm = MPI.COMM_WORLD
        ranks = MPI.Comm_rank(comm)
        #To avoid multiple printing of the same line in parallel
        if ranks != 0
          redirect_stderr(devnull)
          redirect_stdout(devnull)
        end
      end
    end

    @testset "Cases Tests $(typeof(backend))" begin
        airfoil_test(with_debug)
        cylinder_test(with_debug)
        liddriven_test(with_debug)
        taylorgreen_test(with_debug)
    end

end

end #end_module