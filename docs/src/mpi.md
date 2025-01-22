# MPI Run
The code can be run in Message Passing Interface (MPI). 
The code is made in such a way that it can run:
 - in the `REPL`, selecting ` SegregatedVMSSolver.solve(simcase,with_debug)`. It is useful to debug the code and visualize errors.
 - in `MPI`, selecting ` SegregatedVMSSolver.solve(simcase,with_mpi)`

 When running in MPI the code cannot be easily executed in the `REPL`. Instead, one has to run them from a terminal using the [`mpiexecjl`](https://juliaparallel.org/MPI.jl/stable/configuration/#Julia-wrapper-for-mpiexec) script as provided by [`MPI.jl`](https://github.com/JuliaParallel/MPI.jl). 

For example with the command: 

 `mpiexecjl -n 4 julia --project=. run_simulation.jl`

It is possible to launch the tests changing the backend.

## Run in a Cluster
It is possible to run a simulation in a cluster, creating a suitable bash file, example the following `run_sim.sh` file can be used.
The user has to specify where `julia` is installed and to activate a suitable project where there is `SegregatedVMSSolver` in the dependancies. The `run_simulation.jl` is the `julia` file of the simulation itself. Be sure that the number of processors specified in simulation is corresponding to the number of processors used for running the simulation in the cluster.

```bash
#!/bin/sh
export PATH=$HOME/julia-1.11.1/bin/:$PATH

srun julia --project=../../../ -O3 --check-bounds=no -L run_simulation.jl
```

This is a command example to launch the simulation on a single node (`-N1`) a cluster using `slurm`, using 80 cores `-n80`, specifing the partition of 80 cores `-p80CORE` and a total time of 400 hours `-t 400:0:0`.
```bash
sbatch -t 400:0:0 -N 1 -n80 -p80CORE run_sim.sh
```