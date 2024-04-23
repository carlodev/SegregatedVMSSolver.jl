module MPITests

using Test
using MPI
using PartitionedArrays

#project directory
pdir = joinpath(@__DIR__,"..","..",".")

procs = 4
function run_mpi_driver(procs,file)
    mpiexec() do cmd
      
         run(`$cmd -n $procs $(Base.julia_cmd()) --project=$pdir $file`)
      
      @test true
    
  end
end


@testset "Commons MPI" begin
  include(joinpath("..", "CommonsTests", "CommonsTests.jl"))
  run_mpi_driver(procs, CommonsTests.tests_common(with_mpi))
end

@testset "Cases Tests MPI" begin
  include(joinpath("..", "TestsCases", "TestsCases.jl"))
  run_mpi_driver(procs, TestsCases.tests_cases(with_mpi))
end




end