module MPITests

using Test
using MPI

#project directory
pdir = joinpath(@__DIR__,"..","..",".")

procs = 4
function run_driver(procs,file)
    mpiexec() do cmd
      
         run(`$cmd -n $procs $(Base.julia_cmd()) --project=$pdir $file`)
      
      @test true
    
  end
end

run_common = joinpath(@__DIR__,"common_driver_mpi.jl")
run_case = joinpath(@__DIR__,"case_driver_mpi.jl")

run_driver(procs, run_common)
run_driver(procs, run_case)


end