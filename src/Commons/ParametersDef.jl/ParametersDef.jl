module ParametersDef
using Parameters

export SimulationCase
export MeshFileCase
export CartesianCase
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
export ClusterParameters
export ExportParameters
export RestartParameters
 
include("Cases.jl")
include("Params.jl")

end