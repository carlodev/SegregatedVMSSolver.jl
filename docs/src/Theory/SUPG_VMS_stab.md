# SUPG and VMS Stabilization



## Stabilized Equations
Recalling the Galerkin Fromulation of the Navier-Stokes Equations:

``B^G  =   \int_\Omega \left ( R_m \cdot v \right )\;d\Omega + \int_\Omega \left ( R_c \cdot q \right )\;d\Omega``


Using same-order interpolation elements for velocity and pressure (e.g. Q1/Q1)

The SUPG (Streamline Upwind Petrov Galerkin) stabilzation term are:

``B^{SUPG}  =     \int_\Omega\bigg (\tau_m(u \cdot\nabla v +\nabla q)\cdot R_m \bigg )d\Omega+ \int_\Omega \tau_c (\nabla \cdot v)R_c \;d\Omega``

The extra terms added by the VMS:

``B^{VMS1}  = \int_\Omega (u\cdot\nabla v')\cdot(\tau_m R_m) \;d\Omega``

``B^{VMS2}  = -\int_\Omega \bigg (\nabla v\cdot(\tau_m R_m \otimes \tau_m R_m) \bigg )\;d\Omega``


## Stabilization Parameters
The stabilization parameters `\tau_m` and `\tau_c` can be computed accordingly to a `ScalarFormulation()` or a `TensorFormulation()` and it can be specified in the code.

### ScalarFormulation
Typically used with the SUPG. 

``\tau_m = \bigg(\dfrac{2|u|}{h_e} +\dfrac{4\nu}{h_e^2} +\dfrac{2}{dt} \bigg)^{-1}``


``\tau_c = (u\cdot u)\tau_m``

Where `h` is the element dimension. It is the square root of the area for 2D case, and the cubic root of the volume for 3D case.

### TensorFormulation
Typically used with the VMS. 

``\tau_m =\bigg( \dfrac{4}{\Delta t^2} + u\cdot GGu + C_I \nu^2 G:G \bigg)^{-1/2}``

``\tau_c = (\tau_c g\cdot g)^{-1}``

Where G is the inverse of the gradient of the map cell.