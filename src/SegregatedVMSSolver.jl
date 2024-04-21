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



include(joinpath("Commons","SolverOptions.jl"))
include(joinpath("Commons","Interfaces.jl"))

include(joinpath("Commons","ParametersDef","ParametersDef.jl"))
include(joinpath("Commons","VectorsOperations.jl"))

include(joinpath("Commons","AddNewTags.jl"))

include(joinpath("Commons","ModelCreation.jl"))
include(joinpath("Commons","BoundaryConditions.jl"))
include(joinpath("Commons","SpaceConditions.jl"))
include(joinpath("Commons","Equations","Equations.jl"))

include(joinpath("Commons","Restart.jl"))
include(joinpath("Commons","ExportUtility.jl"))

include(joinpath("Commons","InitialConditions.jl"))
include(joinpath("Commons","MatrixCreation.jl"))


include(joinpath("Commons","SolveProblem.jl"))

include("Main.jl")


# include(joinpath("utils","ReadAirfoilResults.jl"))


# include(joinpath("utils","WallDistance.jl"))

include("Exports.jl")

end
