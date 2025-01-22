module ParametersDef

using Parameters
using Gridap
using SyntheticEddyMethod
using LinearAlgebra
using Random

using SegregatedVMSSolver.Interfaces
using SegregatedVMSSolver.SolverOptions

export SimulationCase
export SimulationParameters
export Airfoil
export WindTunnel
export Cylinder
export TaylorGreen
export LidDriven
export VelocityBoundaryCase

export TGVBoundaryConditionsType
export TaylorGreen_Periodic_Parameters
export TaylorGreen_Natural_Parameters
export Natural
export Periodic



export create_new_case
export printstructure
export search_field
export @sunpack
export compute_fluctuation

export UserParameters

export TimeParameters
export PhysicalParameters
export TurbulenceParameters
export SolverParameters
export MeshInfo
export MeshParameters
export TurbulenceDomain
export Internal
export Inlet
export GmshMeshParams
export CartesianMeshParams
export ExportParameters
export InitialParameters


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