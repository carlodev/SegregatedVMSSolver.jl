module WallDistanceTests

using Test
using SegregatedVMSSolver
using SegregatedVMSSolver.WallDistance

mesh_file = "models/DU89_2D_A1_M.msh"
D = 2
u_in = 1.0
Re = 500e3
chord = 1.0
walltag = ["airfoil","wake"]

u_start = get_wall_distance(mesh_file, D, u_in, Re, chord, walltag)

@test u_start <: AbstractVector

end #end module