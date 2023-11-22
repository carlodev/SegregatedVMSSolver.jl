# Boundary Layer Initialization
In order to avoid instabilities and help the convergence of the numerical system this tecnique is adopted. It is based on the resolution of the p-poisson, [Bakker2018](@cite).
The algorithm used to solve the non-linear p-Poisson, equation \eqref{equ:ppoisson} equation resembles the Picard method.

``\nabla \cdot (|\nabla u_p|^{p-2} \nabla u_p) = -1 , x\in \Omega ``
``u_p = 0 , x\in \Omega_D ``

The idea is to compute the wall-normal distance, from the airfoil, for each point of the domain. Then setting a threshold based on the Reynolds number, the velocity in the area close to the airfoil follows a cubic function.

![Velocity initialization](https://carlodev.github.io/SegregatedVMSSolver.jl/dev/Uinit.png)


## How to use it
It is possible to initialize the airfoil simulation using the `WallDistance` module. It works only in serial and for 2D meshes but it is possible to initilize a 3D solution from a 2D results.

```julia
using SegregatedVMSSolver
using SegregatedVMSSolver.WallDistance

mesh_file = "models/DU89_2D_A1_M.msh"
D = 2
u_in = 1.0
Re = 500e3
chord = 1.0
walltag = ["airfoil","wake"]

get_initial_conditions(mesh_file, D, u_in, Re, chord, walltag)
```
This will create a `Initial_Conditions.vtu` which can be open in ParaView and `BoundaryLayerInit.csv` which can be used to restart the simulation.
Only the velocity is initilized, not the pressure.