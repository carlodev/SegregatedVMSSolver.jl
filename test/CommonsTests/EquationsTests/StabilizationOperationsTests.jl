module StabilizationOperationTests
using SegregatedVMSSolver
using Gridap
using Gridap.Fields
using GridapDistributed
using PartitionedArrays
using SegregatedVMSSolver.Equations: compute_stab_coeff,momentum_stabilization, continuity_stabilization



function main(mcase,params,uh,Ω,D::Int64)

    compute_stab_coeff(mcase,params)
    scalar_stab_coeff=compute_stab_coeff(ScalarFormulation(),Ω,D)
    momentum_stabilization(uh, scalar_stab_coeff,mcase )
    continuity_stabilization(uh, scalar_stab_coeff,mcase )

    tensor_stab_coeff=compute_stab_coeff(TensorFormulation(),Ω,D)
    momentum_stabilization(uh, tensor_stab_coeff,mcase )
    continuity_stabilization(uh, tensor_stab_coeff,mcase )
  return true
end


end #end module