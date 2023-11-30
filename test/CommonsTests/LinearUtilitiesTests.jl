module LinearUtilitiesTests

using Test
using SegregatedVMSSolver
using LinearAlgebra
using PartitionedArrays

function main(distribute)
v = rand(10)

vv = SegregatedVMSSolver.create_ũ_vector(v)

@test vv == [v,v,v,v]

v_out = SegregatedVMSSolver.update_ũ(vv)

@test norm(v_out .- v)< 1e-6

new_v = rand(10)

SegregatedVMSSolver.update_ũ_vector!(vv,new_v)

@test vv == [new_v,v,v,v]

end


end