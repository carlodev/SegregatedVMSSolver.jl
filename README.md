# SegregatedVMSSolver

**Documentation**

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://carlodev.github.io/SegregatedVMSSolver.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://carlodev.github.io/SegregatedVMSSolver.jl/dev/)
[![Build Status](https://github.com/carlodev/SegregatedVMSSolver.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/carlodev/SegregatedVMSSolver.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/carlodev/SegregatedVMSSolver.jl/graph/badge.svg?token=GSJB5o6A0m)](https://codecov.io/gh/carlodev/SegregatedVMSSolver.jl)

**Citation**
[![DOI](https://joss.theoj.org/papers/10.21105/joss.07564/status.svg)](https://joss.theoj.org/papers/10.21105/joss.07564)

*SegregatedVMSSolver.jl for solving incompressible Navier-Stokes using stabilized Finite Element Method, in specific Streamline-Upwind Petrov-Galerkin (SUPG) and Variational MultiScale method (VMS)*

![Julia flow](docs/figs/DU89.gif)

## Introduction
The package solves the incompressible Navier-Stokes equations in the Finite Element Framework using SUPG and VMS method. VMS has been originally introduced by [Hughes2000](@cite). In specific, a linearized and segregated version of the SUPG (following the steps illustrated by [Janssens2014](@cite)) and VMS is solved. 

The methods belong to the Large Eddy Simulation (LES) family. The package can solve the `Taylor Green Vortices`, `Lid Driven Cavity Flow` (only 2D), `Cylinder vortex shedding` and and general `Airfoils` (2D and 3D). 

It works fully in parallel. It is specialized for the resolution of flow over airfoils, testing the capability of detecting the Laminar Separation Bubble. It is equipped with some utility `modules` for reading the output files and creating proper initial conditions.

## Installation
The package is registered, so you can install it as:
```julia
using Pkg
Pkg.add(SegregatedVMSSolver)
```

or from the `REPL` just press `]`.

```example
(@1.8) pkg> add SegregatedVMSSolver
```

You can use the most recent release installing it as:
```julia
using Pkg
Pkg.add(url="https://github.com/carlodev/SegregatedVMSSolver.jl")
```

## Suggested software to install

For a complete and smooth experience is suggested to install the free software [ParaView](https://www.paraview.org/) which allows to graphically visualize the results and open `.vtu` and `.pvtu` files.
For creating mesh and physical boundary conditions is suggested to install the free software [`gmsh`](https://gmsh.info/).

## Features
- Implementation of SUPG and VMS formulation for same-order elements for velocity and pressure
- Solve 3D airfoils geometries, time-dependend, fully parallelized code
- Using custom Meshes created with [`gmsh`](https://gmsh.info/). For airfoils the package `AirfoilGmsh.jl` has been developed for speeding up the process
- Solve 2D and 3D cases
- Possibility of choosing the backend thanks to `PartitionedArrays.jl`. It can be run in the `REPL` for debugging or in `MPI`


## Examples

| ![](docs/figs/TGV.gif)  [Taylor-Green Vortices](https://carlodev.github.io/SegregatedVMSSolver.jl/dev/Cases/taylorgreen/) |  ![](docs/figs/Ldx.png)  [Lid Driven Cavity Flow](https://carlodev.github.io/SegregatedVMSSolver.jl/dev/Cases/liddriven/) |
|:-------------:|:-------------:|

| ![](docs/figs/Cyl3D.gif)  [Cylinder Vortex Shedding](https://carlodev.github.io/SegregatedVMSSolver.jl/dev/Cases/cylinder/)  | ![](docs/figs/DU89.gif) [Airfoil](https://carlodev.github.io/SegregatedVMSSolver.jl/dev/Cases/airfoil/)   | 
|:-------------:|:-------------:|


## Parallelization
For parallelization is used `MPI`, and to solve the sparse and distribute numerical systems we use `PETSc`.
The benchmark case is the 2D taylor Green, the time reported here are intended for each time-step. The order of the elements for this simulation is always 2, and the CFL constant at 0.32.

### Strong Parallelization
Strong scalability evaluates how efficiently a parallel code reduces execution time when the problem size remains fixed, but the number of processing units increases. There is a total of `400` elements on each side, leading to `160000` elements and `1920000` dofs in total.
![MPI-strong](docs/figs/STRONG_TGV.png)


### Weak Parallelization
Weak scalability measures how well a parallel code maintains performance when the problem size is kept constant per processor, and the number of processors increases. On each processor there are `50x50` elements, the number of dofs is kept constant at `30K dfos/procs`.
![MPI-weak](docs/figs/WEAK_TGV.png)



## Packages
It relies on the  [`Gridap`](https://github.com/gridap/Gridap.jl) ecosystem. It is also completely written in Julia and allows parallelization. The [`MPI`](https://github.com/JuliaParallel/MPI.jl) and [`PartititionedArrays`](https://github.com/fverdugo/PartitionedArrays.jl) are also at the basis of the parallelization.



## Contributing
It is a collaborative project open to contributions. You can:
- Open a new issue
- Contact the project administator
- Open a PR with the contribution