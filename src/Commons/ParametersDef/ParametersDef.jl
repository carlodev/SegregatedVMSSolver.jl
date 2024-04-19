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
export get_field

export TimeParameters
export PhysicalParameters
export SolverParameters
export MeshParameters
export CartesianMeshParams
export GmshMeshParams
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

include("../AnalyticalSolution.jl")

include("Params.jl")
include("StabilizationStruct.jl")
include("Cases.jl")

end