# MPI Run
The code can be run in Message Passing Interface (MPI). 
The code is made in such a way that it can run:
 - in the `REPL`, selecting ` SegregatedVMSSolver.main(simcase,with_debug)`
 - in `MPI`, selecting ` SegregatedVMSSolver.main(simcase,with_mpi)`

 When running in MPI the code cannot be easily executed in the `REPL`. Instead, one has to run them from a terminal using the [`mpiexecjl`](https://juliaparallel.org/MPI.jl/stable/configuration/#Julia-wrapper-for-mpiexec) script as provided by [`MPI.jl`](https://github.com/JuliaParallel/MPI.jl). 

For example with the command: 

 `mpiexecjl -n 4 julia --project=. run_simulation.jl`

It is possible to launch the tests changing the backend.