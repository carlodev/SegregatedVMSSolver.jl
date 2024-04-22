module ParametersDef

using Parameters

using SegregatedVMSSolver.SolverOptions

export SimulationCase
export SimulationParameters
export Airfoil
export WindTunnel
export Cylinder
export TaylorGreen
export LidDriven
export VelocityBoundaryCase

export printstructure
export search_field
export @sunpack

export UserParameters

export TimeParameters
export PhysicalParameters
export SolverParameters
export MeshInfo
export GmshMeshParams
export CartesianMeshParams
export ExportParameters
export RestartParameters


export StabilizationMethod
export StabilizationFormulation
export StabilizationParameters
export StabilizedProblem
export ScalarStabilization
export TensorStabilization
export ScalarFormulation
export TensorFormulation
export VMS
export SUPG

include("AnalyticalSolution.jl")

include("Params.jl")
include("StabilizationStruct.jl")
include("Cases.jl")

end