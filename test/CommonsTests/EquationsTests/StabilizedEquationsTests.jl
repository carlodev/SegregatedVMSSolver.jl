module StabilizedEquationsTests
using SegregatedVMSSolver
using Gridap
using Gridap.Fields
using GridapDistributed
using PartitionedArrays
using SegregatedVMSSolver.Equations


function main(u_adv,params,simcase)
    vecmat = segregated_equations(u_adv,params,simcase)
    return true
end

end #end module