# Visualization
In this section is explained how to visualize the results.

## Using Paraview
[ParaView](https://www.paraview.org/) allows to graphically visualize the results and open `.vtu` and `.pvtu` files. There are a lot of embedded and advanced tools.

By default, no `.vtu` are created, to save disk space. Creating `Log/PrintSim.txt` allows to monitor the current state of the simulation creating `.pvtu` files. If you cancel the directory, or re-name the file, the code will stop creating `.pvtu` files, it a way to interact in-real-time with the application running.


You can create a Paraview files collecting in sequence all the  `.vtu` to visualize them in temporal sequence. Use the provided api specifying the folder where your `.vtu` are stored.

```julia
using SegregatedVMSSolver.CreateVtu

create_vtu_file("Results_vtu/")
```

!!! info "numeric info" 
    It is suggested to use it at the beginning of the simulation to check the convergence of the simulation and the boundary conditions.
    It may consume a lot of storage space saving all time steps for 3D simulations. 