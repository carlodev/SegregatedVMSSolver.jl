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
export VMS_activation
export Rm_adv_update

include("StabilizationParameters.jl")
include("StabilizationOperations.jl")
include("StabilizedEquations.jl")


end