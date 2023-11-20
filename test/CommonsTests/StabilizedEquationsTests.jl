module StabilizedEquationsTests

using Test
using SegregatedVMSSolver
using Gridap
using GridapDistributed
using PartitionedArrays



function test_stab(rank_partition,distribute,D)
    ν = 0.001
    dt = 0.1
    θ = 0.5
    order = 1
    Cᵢ=[4,36]

parts  = distribute(LinearIndices((prod(rank_partition),)))
domain = (D==2) ? (0,1,0,1) : (0,1,0,1,0,1)
mesh_partition =  (D==2) ?  (4,4) : (4,4,4)
model = CartesianDiscreteModel(parts,rank_partition,domain,mesh_partition)
Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)
params = Dict(:D=>D, :ν=>ν, :dt=>dt, :θ=>θ, :dΩ=>dΩ, :Ω=>Ω,:Cᵢ=>Cᵢ)

vz((x,y)) = (D==2) ? VectorValue(0.0,0.0) : VectorValue(0.0,0.0,0.0)
va((x,y) )= (D==2) ? VectorValue(1.0,2.0) : VectorValue(1.0,2.0,3.0)

reffe =  ReferenceFE(lagrangian,VectorValue{D,Float64},order)
V = TestFESpace(model,reffe,dirichlet_tags="boundary")
U = TrialFESpace(vz,V)
uadv = interpolate_everywhere(va,U)

supg_f = SegregatedVMSSolver.segregated_equations_SUPG!(uadv, params)
@testset "Equations SUPG" for f in supg_f @test typeof(f) <: Function end


vms_f = SegregatedVMSSolver.segregated_equations_VMS!(uadv, params)
@testset "Equations VMS" for f in vms_f @test typeof(f) <: Function end

end


for D in [2,3]
    rank_partition = (D==2) ?  (2,2) : (2,2,2)
    with_debug() do distribute
        test_stab(rank_partition,distribute,D)
    end
end


end


ones(Int64,3)