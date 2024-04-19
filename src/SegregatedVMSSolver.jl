module SegregatedVMSSolver

using Revise
using Gridap
using GridapDistributed
using GridapGmsh
using GridapPETSc
using LinearAlgebra
using NearestNeighbors

using PartitionedArrays
using MPI
using Parameters
using SparseArrays
using FillArrays

using CSV
using DataFrames

using Statistics
using FFTW
using ScatteredInterpolation

using GridapDistributed: Algebra
using Gridap:FESpaces
# using Gridap.Arrays
# using Gridap.CellData



include(joinpath("Commons","SolversOptions.jl"))
include(joinpath("Commons","ParametersDef","ParametersDef.jl"))

include(joinpath("Commons","AddNewTags.jl"))

include(joinpath("Commons","ModelCreation.jl"))
include(joinpath("Commons","BoundaryConditions.jl"))
include(joinpath("Commons","SpaceConditions.jl"))
include(joinpath("Commons","Equations","Equations.jl"))

include(joinpath("Commons","Restart.jl"))
include(joinpath("Commons","InitialConditions.jl"))


# include(joinpath("Commons","ExportUtility.jl"))

# export creation_fe_spaces
# export create_initial_conditions
# export create_PETSc_setup
# export solve_case
include(joinpath("Commons","CommonsProcedures.jl"))


# export create_ũ_vector
# export update_ũ_vector!
# export update_ũ
# include(joinpath("Commons","LinearUtilities.jl"))



# export allocate_Mat_inv_ML
# export inv_lump_vel_mass!
# export initialize_vectors
# export matrices_and_vectors 
# include(joinpath("Commons","MatrixCreation.jl"))

include("Main.jl")


# include(joinpath("utils","ReadAirfoilResults.jl"))


# include(joinpath("utils","WallDistance.jl"))

include("Exports.jl")

end
