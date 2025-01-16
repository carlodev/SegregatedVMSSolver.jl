# Finite Element Method (FEM) for Incompressible Navier-Stokes Equations

The Finite Element Method (FEM) is a widely-used numerical technique for solving partial differential equations (PDEs) that arise in engineering and physics. For the incompressible Navier-Stokes equations, FEM is particularly effective in handling complex geometries and boundary conditions. This section introduces the role of test functions in FEM, the formulation of velocity and pressure spaces, and discusses the issue of Q1/Q1 instabilities.

## Introduction to Test Functions

In FEM, the governing equations are reformulated in their weak form. This involves multiplying the equations by test functions and integrating over the computational domain. For the incompressible Navier-Stokes equations, let:

- ``v``: Test function for velocity, ``v \in V``, where ``V`` is the velocity test function space.
- ``q``: Test function for pressure,  ``q \in Q``, where  ``Q `` is the pressure test function space.

The weak form of the equations is:

``\int_\Omega  (  \frac{\partial u}{\partial t}\cdot v + u \cdot\nabla (u) \cdot v +\nu \nabla (u) \nabla (v) - v\cdot f )\;d\Omega + \int_\Omega  ( \nabla\cdot (u) \cdot q  )\;d\Omega`` 

## Elements and Instabilities

The Q1/Q1 element is a finite element pair where both velocity and pressure are approximated using bilinear shape functions (\( Q1 \) functions) on quadrilateral elements. While computationally simple and efficient, this choice suffers from numerical instabilities due to a violation of the Ladyzhenskaya-Babu\v{s}ka-Brezzi (LBB) or inf-sup condition.

The inf-sup condition ensures a stable coupling between velocity and pressure spaces. Q1/Q1 elements typically fail to satisfy this condition, leading to spurious pressure modes and oscillations in the solution. These instabilities manifest as:

- Unphysical pressure fields with checkerboard patterns.
- Poor convergence of the numerical solution.

### Remedies for Q1/Q1 Instabilities

**Stabilization Techniques**: Introduce stabilization methods like the Streamline Upwind Petrov-Galerkin (SUPG) method to regularize the pressure field.


