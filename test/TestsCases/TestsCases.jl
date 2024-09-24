module TestsCases

using PartitionedArrays
using SegregatedVMSSolver
using Test
using MPI
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions

include("TurbulentAirfoilTest.jl")
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
    turbulent_airfoil_test(backend, true)
    turbulent_airfoil_test(backend, false)
    airfoil_test(backend)
    cylinder_test(backend)
    liddriven_test(backend)
    taylorgreen_test(backend)
  end

end

end #end_module