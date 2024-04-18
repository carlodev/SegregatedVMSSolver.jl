"""
  allocate_Mat_inv_ML(Mat_ML::PSparseMatrix) 

It allocates a zero vector where to store the inverse of the lumped matrix
"""
function allocate_Mat_inv_ML(Mat_ML::PSparseMatrix) 
  return pzeros(Mat_ML.row_partition)
end

function allocate_Mat_inv_ML(Mat_ML::SparseMatrixCSC) 
  l = size(Mat_ML)[1]

  return zeros(l)
end



"""
    set_zeros!(fields::DebugArray)

Set zeros as free values for a field
"""
function set_zeros!(fields::DebugArray)
  for a in fields.items
    a.free_values .= 0.0
  end
end


function set_zeros!(fields::MPIArray)
    fields.item.free_values .= 0.0
end


"""
  inv_lump_vel_mass!(Mat_inv_ML::PVector,Mat_ML::PSparseMatrix)

It computes the lumped matrix, takes the inverse of the diagonal elements.
"""
function inv_lump_vel_mass!(Mat_inv_ML::PVector,Mat_ML::PSparseMatrix)
    values = map(Mat_ML.matrix_partition) do val
        N = maximum(rowvals(val))
    
        V = zeros(N)
        vals = nonzeros(val)
        
        j = 1
        for i in rowvals(val)
            V[i] += vals[j]
            j+1
        end
        V = 1 ./V

    end
    Mat_inv_ML .= PVector(values,Mat_ML.row_partition) 
end


function inv_lump_vel_mass!(Mat_inv_ML::Vector, Mat_ML::SparseMatrixCSC)
  inv_ML_vec = 1 ./ sum(Mat_ML, dims=2)[:,1]
      if !isempty(inv_ML_vec[inv_ML_vec.==Inf])
      error("The matrix ML can not be inverted because after lumping zero values are detected")
  end
  
  Mat_inv_ML.=inv_ML_vec
  
end

"""
  initialize_vectors(matrices::Tuple,uh0,ph0)

It initializes vectors where velocity, pressure, acceleration and all the increments will be stored.
"""
function initialize_vectors(matrices::Tuple,uh0,ph0)
  Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App, Mat_ML, Mat_inv_ML, Mat_S, Vec_Au, Vec_Ap = matrices
  vec_pm = GridapDistributed.change_ghost(get_free_dof_values(ph0), Mat_Aup)
  vec_um = GridapDistributed.change_ghost(get_free_dof_values(uh0), Mat_Auu)


  vec_am = pazeros(Mat_ML)
  vec_sum_pm = pazeros(Mat_Aup)
  Δa_star = pazeros(Mat_Apu)
  Δpm1 = pazeros(Mat_S)
  Δa = pazeros(Mat_Tpu)

  b1 = pazeros(Vec_Au)
  b2 = pazeros(Vec_Ap)
  ũ_vector = create_ũ_vector(vec_um)

  return vec_pm,vec_um,vec_am,vec_sum_pm,Δa_star,Δpm1,Δa,b1,b2,ũ_vector
end



function initialize_matrices_and_vectors(trials,tests, t::Real, u_adv, params)
  return matrices_and_vectors(trials, tests, t::Real, u_adv, params)
end


"""
  matrices_and_vectors(trials, tests, t::Real, u_adv, params)

It updates matrices and vectors
"""
function matrices_and_vectors(trials, tests, t::Real, u_adv, params)

  Tuu,Tpu,Auu,Aup,Apu,App,ML,S,rhs =  segregated_equations(u_adv,params)

    U,P = trials
    V,Q = tests

    Af_Tuu = AffineFEOperator(Tuu,rhs,U(t),V)
    Af_Tpu = AffineFEOperator(Tpu,rhs,U(t),Q)

    Af_Auu = AffineFEOperator(Auu,rhs,U(t),V)
    Af_Aup = AffineFEOperator(Aup,rhs,P(t),V)
    Af_Apu = AffineFEOperator(Apu,rhs,U(t),Q)
    Af_App = AffineFEOperator(App,rhs,P(t),Q)

    Af_ML = AffineFEOperator(ML,rhs,U(t),V)
    Af_S = AffineFEOperator(S,rhs,P(t),Q)

    Mat_Tuu = get_matrix(Af_Tuu)
    Mat_Tpu = get_matrix(Af_Tpu)

    Mat_Auu = get_matrix(Af_Auu)
    Mat_Aup = get_matrix(Af_Aup)
    Mat_Apu = get_matrix(Af_Apu)
    Mat_App = get_matrix(Af_App)

    Mat_ML = get_matrix(Af_ML)
    Mat_S = get_matrix(Af_S)
    
    Vec_Auu = get_vector(Af_Auu)
    Vec_Aup = get_vector(Af_Aup)
    Vec_Apu = get_vector(Af_Apu)
    Vec_App = get_vector(Af_App)

    Mat_inv_ML = allocate_Mat_inv_ML(Mat_ML)
    inv_lump_vel_mass!(Mat_inv_ML,Mat_ML)

    Vec_Ap = Vec_Apu + Vec_App
    Vec_Au = Vec_Auu + Vec_Aup
    return  Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App, Mat_ML, Mat_inv_ML, Mat_S, Vec_Au, Vec_Ap

end
