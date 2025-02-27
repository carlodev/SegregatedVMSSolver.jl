module CommonsTests

using SegregatedVMSSolver
using Test
using PartitionedArrays
using MPI
using Gridap,GridapDistributed


include(joinpath("..","case_test.jl")) 
include(joinpath("CreateProblemTests", "CreateProblemTests.jl"))
include(joinpath("EquationsTests", "EquationsTests.jl"))

include("MatrixCreationTests.jl")
include( "TurbParamsTests.jl")






function tests_common(backend)
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
  
  mcase = create_simulation_test(TaylorGreen, 2)
  params = SegregatedVMSSolver.setup_case(mcase,distribute)
  uh = interpolate_everywhere(VectorValue(0.0,0.0), params[:U])
  ph = interpolate_everywhere(0.0, params[:P](0.0))

  CreateProblemTests.main(distribute)
  EquationsTests.main(uh, params, mcase)
  MatrixCreationTests.main(uh,ph, params, mcase)
  TurbParametersTests.main(distribute)
end #end do distribute

end


tests_common(with_debug)




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