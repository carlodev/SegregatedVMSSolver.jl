# Restart
The package allows to re-start a simulation from a previous stored result by setting:
`:restart=>true`
`mesh_file=>"meshfile.csv"`

It can read a `.csv` file with the headers: `x y z uh_0 uh_1 uh_2 ph`. The `ph` is optional. For each node coordinates `x y z` the velocity in the 3 directions is specified.

It is also possible to use a 2D solution to start a 3D simulation.

