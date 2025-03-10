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
using Gridap.Fields
using GridapDistributed: Algebra
using SegregatedVMSSolver.ParametersDef

export segregated_equations
export coupled_equations

include("StabilizationParameters.jl")
include("StabilizationOperations.jl")
include("StabilizedEquations.jl")
include("CoupledEquations.jl")


end