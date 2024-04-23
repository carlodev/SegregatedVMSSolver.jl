# Simulation Parameters

The simualtion parameters are used to solve the incompressible Navier-Stokes problem.

The Simulation Parameters are defined in `src/Commons/ParametersDef/StabilizationStruct.jl`, where the default values are also provided.

```julia
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef

sprob = StabilizedProblem()
```

The user is creating a `StabilizedProblem` using the default options.
It is possible to specify:
-   method: `VMS()` or `SUPG` where also the order of elements can be specified
-   coeff_method: `TensorFormulation()`, `ScalarFormulation`. Stabilization coefficients `r` and `Ci` can be provided.
-   skew: enable a skew-symmetric formulation

```julia
sprob = StabilizedProblem(method=SUPG(1), coeff_method=TensorFormulation(Ci=[2,16], r=2), skew=true)
```
In this example the `StabilizedProblem` is using a `SUPG` resolution, using first order elements, a method to compute the coeffients based on the tensor formulation, with the supplied coefficients and a skew symmetric formulation of the conservation equations.

