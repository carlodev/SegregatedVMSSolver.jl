# Post Processing

In this section is explained how to visualize the results and use the integrated post processing api for studying airfoils. All the results are saved in the folder `/Results`. 

## Using Paraview
[ParaView](https://www.paraview.org/) which allows to graphically visualize the results and open `.vtu` and `.pvtu` files. There are a lot of embedded and advanced tools.

You can create a Paraview files collecting in sequence all the  `.vtu` to visualize them in temporal sequence. Use the provided api specifing the folder where your `.vtu` are stored.
´´´julia
using SegregatedVMSSolver.CreateVtu

create_vtu_file("Results_vtu/")
´´´

!!! info "numeric info" 
    Creating `Log/PrintSim.txt` allows to monitor the current state of the simulation creating `.pvtu` files.
    It is suggested to use it at the beginning of the simulation to check the convergence of the simulation and the boundary conditions.
    It may consume a lot of storage space saving all time steps for 3D simulations. 

## Using Integrated API
A specific module `ReadAirfoilResults` has been develop. While running an `Airfoil` simulation the code automatically saves the results - for pressure and velocity normal gradient - just for the nodes that are part of the `airfoil` boundary. 

Use the `ReadAirfoilResults` module
```julia
using SegregatedVMSSolver
using SegregatedVMSSolver.ReadAirfoilResults
```

Specify the path where the `.csv` files with the results are stored. 
Specify also the geometrical and physical parameters
```julia
res_path = "Results/"
Re = 500_000
u0 = 1.0
c = 1.0
rho = 1.0
μ = u0*c*rho/Re
α = 1.0
```

Get the nodes and normals
```julia
nodes, normals = get_geometry_info(res_path;α=α)
```
You can get the everage - in time and spanwise direction - of velocity and friction field. It is possible to specify the number of time-step that you want to skip at the beginning of the averaging using the keyword `offset`. It allows to avoid averaging also the initils time-steps where the solution is still evolving.

```julia
Ph = time_space_average_field(res_path, "ph", nodes)
Velocity = time_space_average_field(res_path, "uh", nodes)
Friction = time_space_average_field(res_path, "friction", nodes)
```


Writing the parameters of the simulation allows the code to get the local values for Cp and Cf distiguishing between top and bottom side. 
```julia
cp_top, cp_bottom = extract_Cp(nodes, Ph; u0=u0, rho=rho)
friction_top, friction_bottom = extract_Cf(nodes, Friction, μ; u0=u0, rho=rho)

CL,CD = compute_CL_CD(nodes, cp_top, cp_bottom,
friction_top, friction_bottom; chord=1.0)
```

It is possible to plot the results.
For example

```julia
using Plots
plot(nodes.top.x, cp_top, label = "VMS", color = :red, linewidth = 1.5)
plot!(nodes.bottom.x,cp_bottom, label = false, color = :red, linewidth = 1.5 )
```


# Advanced Post Processing
There are also more powerful API to extract velocity fluctuations, TKE and spectra.

## Time Averaging
In this example we compute the time-average of the velocity field in the volume `"topairfoil"`. The user while creating the mesh has to add this tag to the volume to analyze

```julia
tagname="topairfoil"

topairfoil_nodes = get_nodes(res_path; tagname=tagname)

Vel_avg3D = compute_time_average(res_path; tagname=tagname)
Vel_avg2D = compute_time_span_average(res_path,Vel_avg3D; tagname=tagname)
```

It is possible to extract one of the components of the average velocity
```julia
ux_avg = Vel_avg2D[:,1]
uy_avg = Vel_avg2D[:,2]
uz_avg = Vel_avg2D[:,3]
```
And to plot it on the plane `0.1`

```julia
xgrid,ygrid,velocity_dense = compute_scatter_interp(res_path, ux_avg, 0.1)

using Plots

#it can take a while
scatter(xgrid,ygrid .+ 0.0028;
zcolor=velocity_dense,markerstrokewidth=0, label="U_x average")
```



## Fluctuations

It is possible to create a proble aligned with the `z` axis in the point `xc,yc` and get the fluctuations of that line.

```julia
xc,yc= 0.5,0.1 #it has to be a point of the domain
zprobe, Vel = read_fluctuations(res_path, [xc,yc]; offset=1_000,offend=5_000)
```

It is also possible to compute the Power Spectral Density (PSD)
```julia
dt = 0.001
nt=size(Vel)[1]-1 #Number of time-steps
tn = collect(0.0:dt:dt*nt)
PSD,freqs = compute_PSD(Vel, tn)
```

Computing the TKE in the plane `0.1`
```julia
tke_xy = compute_plane_tke(res_path; zp =0.1)
```



!!! info "OutOfMemory() Error" 
    Processing Large Amount of data can saturate all the RAM available on your system. Try changing machines or reducing the admount of data to process.