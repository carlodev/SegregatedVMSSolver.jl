module MatrixCreation
using Gridap
using Gridap.FESpaces
using GridapDistributed
using SparseArrays
using PartitionedArrays
using LinearAlgebra
using Parameters

using Gridap.Arrays

using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.Equations
using SegregatedVMSSolver.VectorsOperations

export initialize_vectors
export initialize_matrices
export compute_matrices
export update_all_matrices_vectors!


"""
    allocate_Mat_inv_ML(Mat_ML::PSparseMatrix)

Allocate a zero PVector where the inverse of the lumped mass matrix is stored.
"""
function allocate_Mat_inv_ML(Mat_ML::PSparseMatrix)
    return pzeros(Mat_ML.row_partition)
end


"""
    inv_lump_vel_mass!(Mat_inv_ML::PVector, Mat_ML::PSparseMatrix)

Compute the row-sum (lumped) approximation of `Mat_ML`, then store its
reciprocal in `Mat_inv_ML`. Operates locally on each rank.
"""
function inv_lump_vel_mass!(Mat_inv_ML::PVector, Mat_ML::PSparseMatrix)
    values = map(Mat_ML.matrix_partition) do val
        N    = maximum(rowvals(val))
        V    = zeros(N)
        vals = nonzeros(val)
        j    = 1
        for i in rowvals(val)
            V[i] += vals[j]
            j   += 1
        end
        @. V = 1 / V
        V
    end
    Mat_inv_ML .= PVector(values, Mat_ML.row_partition)
end


"""
    initialize_vectors(matrices::Tuple, uh0, ph0)

Allocate the persistent vectors used by the segregated solver: velocity,
pressure, acceleration, their increments, and the RHS buffers `b1` and `b2`.
"""
function initialize_vectors(matrices::Tuple, uh0, ph0)
    Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App,
    Mat_ML, Mat_inv_ML, Mat_S,
    Vec_Auu, Vec_Aup, Vec_Apu, Vec_App, Vec_Au, Vec_Ap = matrices

    vec_pm = GridapDistributed.change_ghost(get_free_dof_values(ph0), Mat_Aup)
    vec_um = GridapDistributed.change_ghost(get_free_dof_values(uh0), Mat_Auu)

    vec_am     = pazeros(Mat_ML)
    vec_sum_pm = pazeros(Mat_Aup)
    Δa_star    = pazeros(Mat_Apu)
    Δpm1       = pazeros(Mat_S)
    Δa         = pazeros(Mat_Tpu)

    b1 = pazeros(Vec_Au)
    b2 = pazeros(Vec_Ap)
    ũ_vector = create_ũ_vector(vec_um)

    return vec_pm, vec_um, vec_am, vec_sum_pm, Δa_star, Δpm1, Δa, b1, b2, ũ_vector
end


function initialize_matrices(u_adv, params, simcase)
    @info "allocating matrices and vectors"
    matrices = allocate_all_matrices_vectors(u_adv, params, simcase)
    @info "matrices and vectors allocated"

    @info "filling matrices and vectors with values"
    @time update_all_matrices_vectors!(matrices, u_adv, params, simcase)
    @info "matrices and vectors updated"

    return matrices
end


function allocate_all_matrices_vectors(u_adv, params, simcase)
    Tuu, Tpu, Auu, Aup, Apu, App, ML, S, rhs = segregated_equations(u_adv, params, simcase)

    @unpack Utn1, Ptn1, tests = params
    V, Q = tests

    Mat_Tuu = allocate_matrix(Tuu, rhs, Utn1, V)
    Mat_Tpu = allocate_matrix(Tpu, rhs, Utn1, Q)

    Mat_Auu, Vec_Auu = allocate_matrix_and_vector(Auu, rhs, Utn1, V)
    Mat_Aup, Vec_Aup = allocate_matrix_and_vector(Aup, rhs, Ptn1, V)
    Mat_Apu, Vec_Apu = allocate_matrix_and_vector(Apu, rhs, Utn1, Q)
    Mat_App, Vec_App = allocate_matrix_and_vector(App, rhs, Ptn1, Q)

    Mat_ML = allocate_matrix(ML, rhs, Utn1, V)
    Mat_S  = allocate_matrix(S,  rhs, Ptn1, Q)

    Mat_inv_ML = allocate_Mat_inv_ML(Mat_ML)
    Vec_Ap     = Vec_Apu + Vec_App
    Vec_Au     = Vec_Auu + Vec_Aup

    return Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App,
           Mat_ML, Mat_inv_ML, Mat_S,
           Vec_Auu, Vec_Aup, Vec_Apu, Vec_App, Vec_Au, Vec_Ap
end


function allocate_matrix(a::Function, rhs::Function, U, V)
    feop = AffineFEOperator(a, rhs, U, V)
    return get_matrix(feop)
end

function allocate_matrix_and_vector(a::Function, rhs::Function, U, V)
    feop = AffineFEOperator(a, rhs, U, V)
    return get_matrix(feop), get_vector(feop)
end


"""
    update_all_matrices_vectors!(matrices, u_adv, params, simcase)

Reassemble all bilinear-form matrices and their associated Dirichlet
contribution vectors in place, then refresh the lumped inverse `Mat_inv_ML`
and the combined Vec_Au, Vec_Ap.
"""
function update_all_matrices_vectors!(matrices::Tuple, u_adv, params, simcase)
    Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App,
    Mat_ML, Mat_inv_ML, Mat_S,
    Vec_Auu, Vec_Aup, Vec_Apu, Vec_App, Vec_Au, Vec_Ap = matrices

    @unpack Utn1, Ptn1, tests = params
    V, Q = tests

    Tuu, Tpu, Auu, Aup, Apu, App, ML, S, _ = segregated_equations(u_adv, params, simcase)

    update_matrix!(Tuu, Mat_Tuu, Utn1, V)
    update_matrix!(Tpu, Mat_Tpu, Utn1, Q)

    update_matrix_vector!(Auu, Mat_Auu, Vec_Auu, Utn1, V)
    update_matrix_vector!(Aup, Mat_Aup, Vec_Aup, Ptn1, V)
    update_matrix_vector!(Apu, Mat_Apu, Vec_Apu, Utn1, Q)
    update_matrix_vector!(App, Mat_App, Vec_App, Ptn1, Q)

    update_matrix!(ML, Mat_ML, Utn1, V)
    update_matrix!(S,  Mat_S,  Ptn1, Q)

    inv_lump_vel_mass!(Mat_inv_ML, Mat_ML)

    # In-place combination, no temporary PVector created.
    @. Vec_Ap = Vec_Apu + Vec_App
    @. Vec_Au = Vec_Auu + Vec_Aup
end


function update_matrix_vector!(a::Function, A::AbstractMatrix, b::AbstractVector, U, V)
    dv = get_fe_basis(V)
    du = get_trial_fe_basis(U)

    mat_contribs = a(du, dv)

    uhd  = zero(U)
    data = collect_cell_matrix_and_vector(U, V, mat_contribs, 0.0, uhd)

    assembler = SparseMatrixAssembler(U, V)
    assemble_matrix_and_vector!(A, b, assembler, data)
end


function update_matrix!(a::Function, A::AbstractMatrix, U, V)
    dv = get_fe_basis(V)
    du = get_trial_fe_basis(U)

    mat_contribs = a(du, dv)
    data         = collect_cell_matrix(U, V, mat_contribs)

    assembler = SparseMatrixAssembler(U, V)
    assemble_matrix!(A, assembler, data)
end


end # module
