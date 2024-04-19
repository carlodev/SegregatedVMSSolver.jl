module Equations
using Gridap
using GridapDistributed
using Parameters
using FillArrays
using LinearAlgebra
using Gridap.CellData

using SegregatedVMSSolver.ParametersDef

export segregated_equations


include("StabilizationParameters.jl")
include("StabilizationOperations.jl")
include("StabilizedEquations.jl")


end