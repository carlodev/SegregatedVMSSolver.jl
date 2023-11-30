module CommonsTests
using SegregatedVMSSolver
using Test
using PartitionedArrays
using MPI
include("InitializeParamsTests.jl")
include("AddNewTagsTests.jl")
include("StabParamsTests.jl")
include("LinearUtilitiesTests.jl")
include("StabilizedEquationsTests.jl")
include("MatrixCreationTests.jl")
include("RestartTests.jl")

function test_common_debug()
  with_debug() do distribute
  InitializeParamsTests.main(distribute)
  AddNewTagsTests.main(distribute)
  StabParamsTests.main(distribute)
  LinearUtilitiesTests.main(distribute)
  StabilizedEquationsTests.main(distribute)
  MatrixCreationTests.main(distribute)
  RestartTests.main(distribute)
  end
end



#with_debug
@testset "Commons Debug" begin
  test_common_debug()
end



pdir = joinpath(@__DIR__,"..","..",".")

procs = 4
function run_driver(procs,file)
    mpiexec() do cmd
      
         run(`$cmd -n $procs $(Base.julia_cmd()) --project=$pdir $file`)
      
      @test true
    
  end
end
  
run_driver(procs, "test/TestMPI.jl")


end