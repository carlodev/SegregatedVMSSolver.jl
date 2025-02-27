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
include("TGV_Natural.jl")
include("TGV_Periodic.jl")


function setup_test_case(simcase, distribute)

  params = SegregatedVMSSolver.setup_case(simcase, distribute)
  @test typeof(params)<: Dict

end


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


    @testset "Create Cases" begin
      setup_test_case(turbulent_airfoil_test(false), distribute)
      setup_test_case(airfoil_test(), distribute)
      setup_test_case(cylinder_test(), distribute)
      setup_test_case(liddriven_test(), distribute)
      setup_test_case(TGV_Periodic_test(2), distribute)
      setup_test_case(TGV_Natural_test(2), distribute)
      if backend == with_debug
        setup_test_case(TGV_Periodic_test(3), distribute)
        setup_test_case(TGV_Natural_test(3), distribute)
      end
  

    end




  end #end backend() do 

  @testset "Solve Cases" begin
    SegregatedVMSSolver.solve(TGV_Periodic_test(2), backend)
    SegregatedVMSSolver.solve(airfoil_test(), backend)

  end


end



end #end_module

