module Equations
using Gridap
using GridapDistributed
using Parameters


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
export segregated_equations

# export compute_stab_coeff
# export momentum_stabilization
# export continuity_stabilization

include("EquationsOperations.jl")
include("StabilizationParameters.jl")
include("StabilizedEquations.jl")


end