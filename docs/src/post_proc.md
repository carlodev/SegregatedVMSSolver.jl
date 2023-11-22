# Post Processing

In this section is explained how to visualize the results and use the integrated post processing api for studying airfoils. All the results are saved in the folder `/Results`. 

## Using Paraview
[ParaView](https://www.paraview.org/) which allows to graphically visualize the results and open `.vtu` and `.pvtu` files. There are a lot of embedded and avanced tools.

!!! info "numeric info" 
    Creating `Log/PrintSim.txt` allows to monitor the current state of the simulation creating `.pvtu` files.

## Using Integrated API
A specific module `ReadAirfoilResults` has been develop. While running an `Airfoil` simulation the code automatically saves the results - for pressure and velocity normal gradient - just for the nodes that are part of the `airfoil` boundary. 

Use the `ReadAirfoilResults` module
```julia
using SegregatedVMSSolver
using SegregatedVMSSolver.ReadAirfoilResults
```

Specify the path where the `.csv` file with the results are stored. Get the nodes, normals.
```julia
res_path = "Results/"
nodes = get_nodes(res_path)
normals = get_normals(res_path)
```

You can get the everage - in time and spanwise direction - of velocity and friction field. It is possible to specify the number of time-step that you want to skip at the beginning of the averaging using the keyword `offset`. It allows to avoid averaging also the initils time-steps where the solution is still evolving.

```julia
Ph = average_field(res_path, "ph", nodes)
Friction = average_field(res_path, "friction", nodes)


```

Writing the parameters of the simulation allows the code to get the local values for fricition and pressure and to distiguish between top and bottom side. 
```julia
Re = 500_000
u0 = 1.0
c = 1.0
rho = 1.0
μ = u0*c*rho/Re
α = 1.0

top_nodesx,bottom_nodesx,top_nodesy,bottom_nodesy,cp_top,cp_bottom,
friction_top,friction_bottom,n_Γ_airfoil_top,n_Γ_airfoil_bottom=extract_airfoil_features(nodes, normals, Ph, Friction; u0=u0, μ=μ, rho=rho, α=α, chord=c)

CL,CD = compute_CL_CD(top_nodesx,bottom_nodesx,top_nodesy,bottom_nodesy,cp_top,cp_bottom,
friction_top,friction_bottom; chord = c)
```

It is possible to plot the results.
For example
```julia
using Plots

plot(top_nodesx, friction_top, label = "VMS", color = :red, linewidth = 1.5)
plot!(bottom_nodesx,friction_bottom, label = false, color = :red, linewidth = 1.5 )

```