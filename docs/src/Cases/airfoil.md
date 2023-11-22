# Airfoil
![SD7003 profile at Reynolds 60000](https://carlodev.github.io/SegregatedVMSSolver.jl/dev/sd7003.png)

It is one of the most complex and intersting case. The user has to create a proper mesh in [`gmsh`](https://gmsh.info/) setting the following physical boundaries:
- `inlet` for the inlet
- `outlet` for the outlet
- `airfoil` for the airfoil walls
- `limits` for the top and bottom boundaries

The velocity at the inlet is incresed from `0.0` arriving to the target value `u_in` at `:t_endramp`. This increase the numeric stability. If `:t_endramp` = `:t0` the velocity at the inlet will be immediately `:u_in`. For numeric stability is better to keep `u_in = 1.0`, then fix the Reynolds and so the viscosity will be automatically computed as: `ν = 1/Reynolds`

The pressure is set `0.0` at the `outlet` section. The velocity on the `limits` is set equal to the one at `inlet`.

## Suggested workflow
3D LES are heavy and it is possible to experience divergence issues. It is suggested to use one of the two initilization tecniques: `Velocity ramping` or `Boundary layer initialization`

### Velocity ramping
By setting `t_endramp > t0` automatically the code will create an inlet velocity which will increase linearly in time up to the `u_in` target value.

```julia
uin(t) = (t < t_endramp) ? (u_in .*(t/t_endramp)) : u_in
```


### Boundary Layer initialization
Details are provided in the dedicated section.

