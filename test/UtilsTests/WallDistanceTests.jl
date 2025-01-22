module WallDistanceTests

using Test
using DataFrames
using SegregatedVMSSolver.WallDistance

mesh_file = joinpath(@__DIR__, "../../models", "DU89_2D_A1_M.msh")

D = 2
u_in = 1.0
Re = 500e3
chord = 1.0
walltag = ["airfoil","wake"]

df_start = get_initial_conditions(mesh_file, u_in, Re, walltag; D=D, chord=chord, p=4)
@test typeof(df_start) <: DataFrame


end #end module