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
using Gridap.Arrays
using Gridap.CellData


include("Main.jl")
include(joinpath("ParametersDef","ParametersDef.jl"))

export create_model
include(joinpath("Commons","ModelCreation.jl"))

export create_boundary_conditions
include(joinpath("Commons","BoundaryConditions.jl"))

export creation_fe_spaces
include(joinpath("Commons","SpaceConditions.jl"))


export init_params
export verifykey
include(joinpath("Commons","InitializeParameters.jl"))

include(joinpath("Commons","ExportUtility.jl"))

export print_model
export creation_fe_spaces
export create_initial_conditions
export create_PETSc_setup
export solve_case
include(joinpath("Commons","CommonsProcedures.jl"))

export create_new_tag!
export add_new_tag!
export add_centre_tag!
include(joinpath("Commons","AddNewTags.jl"))


export create_ũ_vector
export update_ũ_vector!
export update_ũ
include(joinpath("Commons","LinearUtilities.jl"))

include(joinpath("Equations","Equations.jl"))


export petsc_options_default
export petsc_options
include(joinpath("Commons","SolversOptions.jl"))

export allocate_Mat_inv_ML
export inv_lump_vel_mass!
export initialize_vectors
export matrices_and_vectors 
include(joinpath("Commons","MatrixCreation.jl"))

export create_search_tree
export restart_uh_field
export restart_ph_field
include(joinpath("Commons","Restart.jl"))




include(joinpath("utils","ReadAirfoilResults.jl"))


include(joinpath("utils","WallDistance.jl"))


#TaylorGreen
include(joinpath("TaylorGreen","TaylorGreen.jl"))

#LidDriven
include(joinpath("LidDriven","LidDriven.jl"))

#Cylinder
include(joinpath("Cylinder","Cylinder.jl"))

#Airfoil
include(joinpath("Airfoil","Airfoil.jl"))

#WindTunnel
include(joinpath("WindTunnel","WindTunnel.jl"))
end
