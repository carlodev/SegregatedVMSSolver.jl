module Equations
using Gridap
using GridapDistributed
using Parameters
using FillArrays
using LinearAlgebra
using PartitionedArrays

using Gridap.CellData
using Gridap.Arrays
using Gridap.CellData
using GridapDistributed: Algebra
using SegregatedVMSSolver.ParametersDef

export segregated_equations


include("StabilizationParameters.jl")
include("StabilizationOperations.jl")
include("StabilizedEquations.jl")


end