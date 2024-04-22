module CreateProblem

using Gridap
using GridapDistributed
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
export create_initial_conditions

include("AddNewTags.jl")
include("Restart.jl")


include("ModelCreation.jl")
include("BoundaryConditions.jl")
include("SpaceConditions.jl")
include("InitialConditions.jl")

end #end module