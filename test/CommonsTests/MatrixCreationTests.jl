module MatrixCreationTests

using Test
using SegregatedVMSSolver
using LinearAlgebra
using PartitionedArrays, SparseArrays

N = 1000
Mat = sprand(N,N,1e-3) + sparse(Matrix(2.0I, N, N))
Mat_inv = SegregatedVMSSolver.allocate_Mat_inv_ML(Mat)
@test Mat_inv == zeros(N)

SegregatedVMSSolver.inv_lump_vel_mass!(Mat_inv, Mat)

b = ones(N)

@test norm(Mat_inv .* b .- Mat\b) ./ N < 1.0

end