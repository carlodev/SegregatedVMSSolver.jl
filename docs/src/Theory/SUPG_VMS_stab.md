# SUPG and VMS Stabilization


Recalling the Galerkin Fromulation of the Navier-Stokes Equations:

``B^G  =   \int_\Omega \left ( R_m \cdot v \right )\;d\Omega + \int_\Omega \left ( R_c \cdot q \right )\;d\Omega``


Using same-order interpolation elements for velocity and pressure (e.g. Q1/Q1)

The SUPG (Streamline Upwind Petrov Galerkin) stabilzation term are:

``B^{SUPG}  =     \int_\Omega\bigg (\tau_m(u \cdot\nabla v +\nabla q)\cdot R_m \bigg )d\Omega+ \int_\Omega \tau_c (\nabla \cdot v)R_c \;d\Omega``

The extra terms added by the VMS:

``B^{VMS1}  = \int_\Omega (u\cdot\nabla v')\cdot(\tau_m R_m) \;d\Omega``

``B^{VMS2}  = -\int_\Omega \bigg (\nabla v\cdot(\tau_m R_m \otimes \tau_m R_m) \bigg )\;d\Omega``


