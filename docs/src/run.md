# Package usage
The package allows the user to set a wide variety of options.
Problem Settings:
- `:N` = number of divisions for each dimension.
- `:D` = dimension. It can be 2 or 3.
- `:order` = order of the elements. At the moment just the order 1 is tested.
- `:case` it can be `"TaylorGreen", "LidDriven", "Cylinder", "Channel", "Airfoil"`.
- `:u_in` the inlet velocity for `"Airfoil"` and `"Cylinder"`, or the lid velocity for `"LidDriven`
- `:c` chord length in the `"Airfoil"` case, or dimension of lid for `"LidDriven"`. It is used to compute the viscosity `:ν` from the Reynolds and velocity
- `:Re` Reynolds number. 
- `:ν` kinematic viscosity. It can be overwritten in order to satisfy the Reynolds, in this case a warning informs the user.
- `:ρ` density. It used just to compute the force. The advice is to keep it `1.0` and just set the Reynolds.

Time settings
- `:t0` starting time.
- `:dt` time step length.
- `:tF` end time.
- `:t_endramp` for high reynolds cases, like airfoils and lid driven, for improving numeric stability the inlet velocity (or the lid velocity) are increased from 0 up the desired value in the time between :t0 and :t_endramp. If :t0 = :t_endramp there is no ramping.

Ode Settings
- ``:θ`` parameter required for time integration.  ``:θ = 0.5`` allows to have 2nd order accuracy on velocity. For the pressure is alwys used a fully implicit method.

Numeric Settings
- `:method` can be `:SUPG` or `:VMS`
- `:Cᵢ` is a vector containing stabilization coefficients used for the `:VMS`. The suggested values are `[4,36]`, [Trofimova2009](@cite)
- `:options` the settings for the petsc solver. The function [`petsc_options`](@ref) can be used.

Print Settings
- `:printmodel` can be `true` or `false`. If `true` mesh is saved as a .pvtu file.
- `:printinitial` can be `true` or `false`. If `true` saves the flowfield at `t0`. It is useful when restarting from a previous solution.
- `:benchmark`  can be `true` or `false`. If `true` it does not print the solution, and it gives the time needed for computing the iteration form the 2nd till the end. The first iteration is not taken into account for computing the time because of precompilation.
- `:name_tags` specify the boundary names to export, eg: `:name_tags=>["airfoil","inlet"]`


Mesh Settings
- `:mesh_file` is a string with the name of the `.msh` that can be read. By default it points to the folder `/models` of the package.



Partitioning Settings
- `:rank_partition` tuple which set the number of division on each axes to split the geometry to the processors.  For non cartesian problems it does not matter how the cores are split into the dimensions as long as `prod(rank_partition)` is equal to the `MPI` ranks.
- `backend` can be `with_debug()` (can be run in the `REPL`) or `with_mpi()`. 

Restart Settings
- `:restart` can be `true` or `false`. If `false` the initial conditions are computed internally using `:u_in` or analytical solution (`"TaylorGreen"`). 
- `:restart_file` is used only if `:restart`is true. It is a `.csv` file created from ParaView using the SpreadSheet. It has the list of fo velocity and pressure in each node. It is better to run `clean grid` in Paraview before for get rid of duplicate points.

Monitor Settings
- `:log_dir` select the folder path where the simulation while running check if there is the file `PrintSim.txt`. If the file exists the simulation will start creating `.pvtu` files at the end of every time step. The contenent of the file does not matter. By default the directory is `/Log`





```Example
using SegregatedVMSSolver
using Parameters,PartitionedArrays

  
   params = Dict(
         :N => 32,
         :D => 2, #Dimension
         :order => 1, 
         :t0 => 0.0,
         :dt => 0.25,
         :tF => 1.0,
         :case => "TaylorGreen",
         :θ => 0.5,
         :u_in=> 1.0,
        
         :backend => with_debug,  #or with_mpi()
         :rank_partition=>(2,2),
         :ν => 0.001,
         :petsc_options => petsc_options_default(),
         :method=>:VMS,
         :Cᵢ => [4, 36],
    
         :t_endramp=>0.0,
         :mesh_file => " ",
         :Re=> 1000,
         
         :c=> 1.0,
 
   )
   
   
 SegregatedVMSSolver.main(params) 

```

The example above run in the `REPL` emulating a parallel run over 4 processors (you can see it by the options `rank_partition`). 


!!! info "numeric info" 
    Changing the backend to `with_mpi()` allows it to run in MPI. `with_debug` is useful to debug the code.

!!! info "numeric info" 
    Creating `Log/PrintSim.txt` allows to monitor the current state of the simulation. 