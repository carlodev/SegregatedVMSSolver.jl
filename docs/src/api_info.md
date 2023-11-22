# Index

## Initialize Parameters
```@docs
init_params
verifykey
```

## Commons Procedures
```@docs
print_model
creation_fe_spaces
create_initial_conditions
create_PETSc_setup
solve_case
```

## AddNewTags
```@docs
create_new_tag!
add_new_tag!
add_centre_tag!
```

## Stabilization Parameters
```@docs
h_param
G_params
compute_d
compute_G
compute_GG
compute_gg
```

## Linear Utilities
```@docs
create_ũ_vector
update_ũ_vector!
update_ũ
```


## Stabilized Equations
```@docs
cconv
segregated_equations_SUPG!
segregated_equations_VMS!
```

## SolversOptions
```@docs
petsc_options
```

## MatrixCreation
```@docs
allocate_Mat_inv_ML
inv_lump_vel_mass!
initialize_vectors
matrices_and_vectors 
```

## Restart
```@docs
find_idx
uh_restart
ph_restart
restart_uh_field
restart_ph_field
```


```@meta
CurrentModule = SegregatedVMSSolver.ExportUtility
```
## SegregatedVMSSolver.ExportUtility

```@autodocs
Modules = [ ExportUtility,]
```



```@meta
CurrentModule = SegregatedVMSSolver.ReadAirfoilResults
```
## SegregatedVMSSolver.ReadAirfoilResults

```@autodocs
Modules = [ ReadAirfoilResults,]
```


```@meta
CurrentModule = SegregatedVMSSolver.WallDistance
```
## SegregatedVMSSolver.WallDistance

```@autodocs
Modules = [ WallDistance,]
```