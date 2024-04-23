# User Parameters

The simulation parameters are written step by step by the user.
The Physical Parameters are defined in `src/Commons/ParametersDef/Params.jl`, where the default values are also provided.

### TimeParameters
```julia
using SegregatedVMSSolver
using PartitionedArrays,MPI
using SegregatedVMSSolver.ParametersDef
t0 = 0.0
dt = 0.1
tF = 1.0

timep = TimeParameters(t0,dt,tF)
```
The user can also specify a `time_window::Tuple{Float64,Float64}` where to compute the average flow-field.
The user can specify a `t_endramp`, the boundary-velocity is then increased from 0.0 up to the target value following a linear law from `t0` up^to `t_endramp`, it increases stability.


### PhysicalParameters
```julia
Re = 1000
c= 1.0
u_in = 1.0
physicalp = PhysicalParameters(Re=Re,c=c,u_in=u_in)
```
The viscosity `ν` is computed automatically in this way

### SolverParameters
```julia
solverp = SolverParameters()
```
In this example we are using the default parameters.
It is possible to supply more detailed options, as reported in the file `src/Commons/SolverOptions.jl`.

```julia
options = petsc_options(vel_ksp="gmres", vel_pc="ilu")
solverp = SolverParameters(petsc_options=options)
```

!!! info "Using HYPRE" 
    For solving the `Airfoil` case, it is suggested to use the option `petsc_options_airfoil()` which requires the installation of `hypre` when installing PETSc


### MeshParameters
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

```julia
exportp = ExportParameters(printinitial=false,printmodel=false,name_tags=["airfoil"], fieldexport=[["uh"]])
```
In this example, the velocity field (or better, the gradient of the velocity field) corresponding to the `airfoil` tag is exported at each time step.


!!! info "Export `.vtu` files" 
    Exporting `.vtu` files can rapidly fill the disk space. By default, a `.vtu` file is exported every 100 time steps. You can enable a continuous exporting of `.vtu` files creating the following file : `Log/PrintSim.txt`


### RestartParameters
```julia
restartp = RestartParameters()
```
The creation of RestartParameters is optional is not enabled by default.
The user can specify a `.csv` file path containing the columns `Points_0,Points_1,Points_2,uh_0,uh_1,uh_2`.
