---
title: 'SegregatedVMSSolver.jl: Linearized and Segregated Stabilized Method for Large Eddy Simulation in Julia'
tags:
  - Julia
  - Turbulence
  - Large Eddy Simulation
  - Variational MultiScale Method
  - VMS
  - Streamline-Upwind/Petrov-Galerkin
  - SUPG
authors:
  - name: Carlo Brunelli
    orcid: 0000-0002-2873-6293
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Bart Janssens
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: Mechanical Engineering Department, Royal Military Academy, Belgium
   index: 1
date: 30 November 2023
bibliography: paper.bib

---

# Summary
Large-Eddy Simulation (LES) is a family of mathematical techniques that perform high-fidelity simulations in Computational Fluid Dynamics (CFD). They can simulate turbulent flows by numerically solving the Navier-Stokes equations. Using filtering operation, LES focuses on the larger length scale, while the effects of small scales (subgrid scales) are modeled. The Variational MultiScale (VMS) and Streamline-Upwind/Petrov-Galerkin (SUPG) methods belong to the family of stabilized methods. The mathematical problem is solved using the Finite Element Method (FEM) framework. They allow to simulate complex flows. The package has been developed with the primary aim of studying the Laminar Separation Bubble (LSB) at low-Reynolds regime on the suction side of the airfoils.

# Statement of need
`SegregatedVMSSolver.jl` is a package that aims to perform high-fidelity simulations of complex flow phenomena in the incompressible regime using the VMS and SUPG method. VMS has been originally introduced by [@Hughes:2000]. The linearization adopted has been proposed in the SUPG method by [@Banyai:2016].

The package is using `Gridap.jl`[@Verdugo:2022],[@Badia:2022] package to implement the mathematical model of FEM. It uses `GridapDistributed.jl`[@BadiaD:2022] and [`PartitionedArrays.jl`](https://github.com/fverdugo/PartitionedArrays.jl) to be able to run on multicore - using MPI protocol. The [`GridapPETSc.jl`](https://github.com/gridap/GridapPETSc.jl) package is used to solve the final linear system. It solves a Linearized and Segregated version of VMS (LS-VMS) and SUPG. It uses the $\theta$ method to solve the time-marching problem. Different 2D test cases are provided: Taylor-Green vortices, lid-driven, and the vortex shedding over a cylinder. The package is actually dedicated to study airfoils aerodynamic features low-Reynolds regime (up to 500 000). It has a suite of tools for post-processing the results, performing time and spanwise averaging. The user can load own airfoil meshes also usign the package [`AirfoilGmsh.jl`](https://github.com/carlodev/AirfoilGmsh.jl).

Different software packages have been developed in the field of fluid dynamics in Julia: [Julia Packages fluid-dynamics](https://juliapackages.com/c/fluid-dynamics) showing a growing interest in the field by the Julia community.

## Results

![Taylor-Green vortices, velocity x direction field](images/TGx.png){ width=50%  }

![lid-driven cavity flow, velocity field, Reynolds 10 000](images/Ldx.png){ width=50%  }

![Cylinder velocity flow-field, Reynolds 1000](images/Cyx.png){ width=50%  }

![Velocity field on DU89 airfoil at Reynolds 250 000, angle of attack 1°](images/DU89U.png){ width=50%  }

Figure \ref{fig:cfsd7003} shows the comparison of the time averaged results obtained using the LS-VMS, compared with the results obtained by Calderer et al. [@Calderer:2013] and Galbraith [@Galbraith:2008].

![Friction coefficient on the suction side of the sd7003 airfoil, Reynolds 60 000, angle of attack 4°,\label{fig:cfsd7003}](images/VMS7003s.pdf){ width=50%  }


# Package Features
- Support 2-dimensional and 3-dimensional geometries
- It solved a time-dependent problem
- It can run in parallel using MPI
- Velocity ramping
- It can capture the laminar-to-turbulent transition
- It has utilities for time and space averaging
- It supports advanced boundary layer initialization
- Real-time simulation interaction
- The solution can be visualized in ParaView
- Restart the simulation from a specific saved time-step

# References
