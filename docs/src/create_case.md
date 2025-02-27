# Create Your Own Case
The code itself allows the user to create custom cases for their needs.

## Box case
```julia
using SyntheticEddyMethod
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using SegregatedVMSSolver.CreateProblem


create_new_case(:Box)
end 
```

The user can create specific boundary conditions for the case. The function requires 3 outputs.

```julia
function SegregatedVMSSolver.CreateProblem.boundary_velocities(simcase::Box)   
    
    
    u_inlet(x,t) = VectorValue(1.0,0.0)
    u_inlet(t::Real) = x -> u_inlet(x,t)

    u_wall(x,t) = VectorValue(0.0,0.0)
    u_wall(t::Real) = x -> u_wall(x,t)

    pressure_outlet(x,t) = 0.0
    pressure_outlet(t::Real) = x -> pressure_outlet(x,t)
    
    return u_inlet,u_outlet,pressure_outlet
end
```

And then, associating the function to the boundaries

```julia
function SegregatedVMSSolver.CreateProblem.create_boundary_conditions(simcase::Box, u_inlet,u_outlet,pressure_outlet) 
        u_diri_tags=["inlet","limits"]
        u_diri_values = [u_inlet,u_outlet]
        p_diri_tags=["outlet"]
        p_diri_values = [pressure_outlet]
        return u_diri_tags,u_diri_values,p_diri_tags,p_diri_values
end 
```


