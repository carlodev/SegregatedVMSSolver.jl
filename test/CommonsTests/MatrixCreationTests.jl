module MatrixCreationTests

using Test
using SegregatedVMSSolver
using LinearAlgebra
using PartitionedArrays, SparseArrays
using Gridap, GridapDistributed
using Gridap:FESpaces

# N = 1000
# Mat = sprand(N,N,1e-3) + sparse(Matrix(2.0I, N, N))
# Mat_inv = SegregatedVMSSolver.allocate_Mat_inv_ML(Mat)
# @test Mat_inv == zeros(N)

# SegregatedVMSSolver.inv_lump_vel_mass!(Mat_inv, Mat)

# b = ones(N)

# @test norm(Mat_inv .* b .- Mat\b) ./ N < 1.0


function test_matrix(rank_partition,distribute,D)
    ν = 0.001
    order = 1


parts  = distribute(LinearIndices((prod(rank_partition),)))
domain = (D==2) ? (0,1,0,1) : (0,1,0,1,0,1)
mesh_partition =  (D==2) ?  (4,4) : (4,4,4)
model = CartesianDiscreteModel(parts,rank_partition,domain,mesh_partition)
Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)

vz((x,y)) = (D==2) ? VectorValue(1.0,0.0) : VectorValue(1.0,0.0,0.0)
va((x,y) )= (D==2) ? VectorValue(1.0,2.0) : VectorValue(1.0,2.0,3.0)

reffe =  ReferenceFE(lagrangian,VectorValue{D,Float64},order)
V = TestFESpace(model,reffe,dirichlet_tags="boundary")
U = TrialFESpace(vz,V)

rhs(v) = 0.0
Au(u, v) = ∫(ν * ∇(v) ⊙ ∇(u) )dΩ 
MAu = get_matrix(AffineFEOperator(Au,rhs,U,V))
@test typeof(MAu) <: PSparseMatrix

Mat_inv = allocate_Mat_inv_ML(MAu)
SegregatedVMSSolver.inv_lump_vel_mass!(Mat_inv, MAu)

@test typeof(Mat_inv)<:PVector

end

function main(distribute)
for D in [2,3]
    rank_partition = (D==2) ?  (2,2) : (2,2,1)
        test_matrix(rank_partition,distribute,D)   
end
end

end