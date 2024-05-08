# Restart
The package allows to re-start a simulation from a previous stored result by setting:
`:restart=>true`
`mesh_file=>"meshfile.csv"`

It can read a `.csv` file with the headers: `Points_0 Points_1 Points_2 uh_0 uh_1 uh_2 ph`. The `ph` is optional. For each node coordinates `x y z` the velocity in the 3 directions is specified. Check in folder `restarts` for examples.

It is also possible to use a 2D solution to start a 3D simulation.

