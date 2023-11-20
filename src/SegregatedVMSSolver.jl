module SegregatedVMSSolver

using Revise
using Gridap
using GridapDistributed
using GridapGmsh
using GridapPETSc
using LinearAlgebra

using PartitionedArrays
using MPI
using Parameters
using SparseArrays
using FillArrays

using CSV
using DataFrames

using GridapDistributed: Algebra
using Gridap:FESpaces
using Gridap.Arrays
using Gridap.CellData

include("Main.jl")
include(joinpath("Commons","Initialize_Parameters.jl"))
include(joinpath("Commons","CommonsProcedures.jl"))
include(joinpath("Commons","AddNewTags.jl"))
include(joinpath("Commons","StabParams.jl"))
include(joinpath("Commons","LinearUtilities.jl"))
include(joinpath("Commons","StabilizedEquations.jl"))
include(joinpath("Commons","SolversOptions.jl"))
include(joinpath("Commons","MatrixCreation.jl"))
include(joinpath("Commons","Restart.jl"))
include(joinpath("Commons","ExportUtilities.jl"))

#TaylorGreen
include(joinpath("TaylorGreen","TaylorGreen.jl"))

#LidDriven
include(joinpath("LidDriven","LidDriven.jl"))

#Cylinder
include(joinpath("Cylinder","Cylinder.jl"))

#Airfoil
include(joinpath("Airfoil","Airfoil.jl"))
end
