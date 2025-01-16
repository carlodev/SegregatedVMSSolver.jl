# The Incompressible Unsteady Navier-Stokes Equations
The package ´SegregatedVMSSolver´ solves the unsteady incompressible Navier-Stokes. The conservation equations takes the following form:

Mass Conservation

``R_c = \nabla\cdot(u)=0 ``

Momentum Conservation

``R_m = \frac{\partial u}{\partial t} + u\cdot\nabla(u) +\nabla(p) -\nu \Delta(u) = f``


Where:
-  ``u`` is the velocity field
-  ``p`` is the pressure (normalized with density, assumed to be unitary)
-  ``\nu`` is the kinematic viscosity of the fluid
- ``f`` is the momentum source term 

