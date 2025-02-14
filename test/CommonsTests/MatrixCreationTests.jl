module MatrixCreationTests

using Test
using SegregatedVMSSolver.ParametersDef
using LinearAlgebra, Parameters
using PartitionedArrays, SparseArrays
using Gridap, GridapDistributed
using Gridap: FESpaces

using SegregatedVMSSolver.MatrixCreation

using SegregatedVMSSolver.MatrixCreation: allocate_all_matrices_vectors


function main(u_adv, p0, params, simcase)

    @testset "Matrix Creation Tests" begin
        @unpack trials, tests = params
        U, P = trials

        @sunpack t0, dt, save_sim_dir = simcase


        Ut0 = U(t0)
        Pt0 = P(t0)

        Ut0_1 = U(t0 + dt)
        Pt0_1 = P(t0 + dt)

        merge!(params, Dict(:Utn => Ut0, :Ptn => Pt0, :Utn1 => Ut0_1, :Ptn1 => Pt0_1))


        matrices = initialize_matrices(u_adv, params, simcase)
        matrices = allocate_all_matrices_vectors(u_adv, params, simcase)
        update_all_matrices_vectors!(matrices, u_adv, params, simcase)
        vec = initialize_vectors(matrices, u_adv, p0)

        @test typeof(matrices) <: Tuple
        @test typeof(vec) <: Tuple
    end

end


end