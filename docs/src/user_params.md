# User Parameters

The simulation parameters are written step by step by the user. In this section a brief explanation of the main parameters is provided.
The Physical Parameters are defined in `src/Commons/ParametersDef/Params.jl`, where the default values are also provided.

### TimeParameters
```julia
using SegregatedVMSSolver
using PartitionedArrays,MPI
using SegregatedVMSSolver.ParametersDef
t0 = 0.0
dt = 0.1
tF = 10.0

timep = TimeParameters(t0,dt,tF)
```
The user can also specify a `time_window::Tuple{Float64,Float64}` where to compute the average flow-field.
The user can specify a `t_endramp`, the boundary-velocity is then increased from 0.0 up to the target value following a linear law from `t0` up to `t_endramp`, it increases stability.

```julia
TimeParameters(t0=t0,dt=dt,tF=tF, t_endramp=2.0, time_window=(5.0,10.0))
```

### PhysicalParameters
```julia
Re = 1000
c= 1.0
u_in = [1.0,0.0]
physicalp = PhysicalParameters(Re=Re,c=c,u_in=u_in)
```
The viscosity `ν` is computed automatically in this way

```julia
physicalp.ν
```

In case of 3D simulation, the initial velocity should be a 3-elements vector.
```julia
u_in = [1.0,0.0,0.0]
```



### SolverParameters
```julia
solverp = SolverParameters()
```
In this example we are using the default parameters.
It is possible to supply more detailed options, as reported in the file `src/Commons/SolverOptions.jl`.

```julia
options = petsc_options(vel_ksp="gmres", vel_pc="gamg")
solverp = SolverParameters(petsc_options=options)
```

The user can specify  `M` the maximum number of internal interations for each time-step and the interval to consider after updating the matrices `matrix_freq_update`.

```julia
solverp = SolverParameters(M=1, matrix_freq_update=10)
```


!!! info "Using HYPRE" 
    For solving the `Airfoil` case, it is suggested to use the option `petsc_options_airfoil()` which requires the installation of `hypre` when installing PETSc


### MeshParameters
For Cartesian cases (Lid-Driven and Taylor-Green)
```julia
D = 2
rank_partition = (2,2)
meshp= MeshParameters(rank_partition,D; N=50,L=0.5)
```
The mesh information are supplied. `D` is the dimension of the problem - 2 or 3. In the example a sqare carthesian mesh is created, with 50 divisions on each side, extending from (-0.5,-0.5) up to (0.5,0.5). The rank_partition is a `Tuple`specifing how to split the domainin the different directions. In this case 4 processors are necessary. 

```julia
airfoil_mesh_file = joinpath(@__DIR__,"..", "..", "models", "DU89_2D_A1_M.msh")
meshp= MeshParameters(rank_partition,D,airfoil_mesh_file)
```
In this example the mesh is created from a `.msh` file created using `Gmsh`.


### ExportParameters
```julia
exportp = ExportParameters()
```
It is using the default values.

The user can choose to export the initial conditions and/or the mesh. It is useful in preliminary stage of simulation design.
```julia
exportp = ExportParameters(printinitial=true,printmodel=true)
```
The package allows to export the value of velocity, pressure or friction from specific targets on the domain.
Example 

```julia
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh"]])
```

!!! info "Export `.vtu` files" 
    Exporting `.vtu` files can rapidly fill the disk space. By default, a `.vtu` file is exported every 100 time steps. You can enable a continuous exporting of `.vtu` files creating the following file : `Log/PrintSim.txt`


### InitialParameters
```julia
intialp = InitialParameters()
```

The user can specify and initial velocity in the whole domain:
```julia
intialp = InitialParameters(u0=[0.5,0.1])
```
The package allows to re-start a simulation from a previous stored result.
It can read a `.csv` file with the headers: `Points_0 Points_1 Points_2 uh_0 uh_1 uh_2 ph`. The `ph` is optional. For each node coordinates `x y z` the velocity in the 3 directions is specified. Check in folder `restarts` for examples of files.
It is also possible to use a 2D solution to start a 3D simulation.

```julia
airfoil_restart_file = joinpath(@__DIR__,"..", "..", "restarts", "BL_DU89_2D_A1_M.csv")
intialp = InitialParameters(airfoil_restart_file)
```

The initial file can also be created using the advanced boundary layer initialization. 