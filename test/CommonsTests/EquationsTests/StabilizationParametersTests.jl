module StabilizationParametersTests
using SegregatedVMSSolver
using Gridap
using Gridap.Fields
using GridapDistributed
using PartitionedArrays
using SegregatedVMSSolver.Equations: h_param, G_params


function main(Ω,D::Int64)
    h_param(Ω, D)
    G_params(Ω, D)

    return true
end

end #end module