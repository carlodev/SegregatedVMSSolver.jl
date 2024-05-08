module MPITests

using Test
using MPI
using PartitionedArrays

#project directory
pdir = joinpath(@__DIR__,"..","..",".")

mpidir = @__DIR__

procs = 4
function run_mpi_driver(procs,file)
    mpiexec() do cmd
      
         run(`$cmd -n $procs $(Base.julia_cmd()) --project=$pdir $(joinpath(mpidir,file))`)
      
      @test true
    
  end
end


run_mpi_driver(procs, "mpi_test.jl")





end