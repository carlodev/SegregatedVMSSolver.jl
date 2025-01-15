module CreateProblem

using Gridap
using GridapDistributed
using GridapGmsh
using Gridap.Arrays
using Gridap.CellData

using Parameters
using PartitionedArrays
using CSV
using DataFrames
using NearestNeighbors

using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.ExportUtility

export create_model
export create_boundary_conditions
export creation_fe_spaces
export creation_∇fe_space
export create_initial_conditions
export create_initial_∇uh
include("AddNewTags.jl")
include("Restart.jl")


include("ModelCreation.jl")
include("BoundaryConditions.jl")
include("SpaceConditions.jl")
include("InitialConditions.jl")

end #end module