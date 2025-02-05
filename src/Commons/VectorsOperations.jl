module VectorsOperations

using Gridap
using GridapDistributed
using PartitionedArrays
using SparseArrays
using Parameters
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef

export create_ũ_vector
export update_ũ_vector!
export update_ũ
export pazeros
export set_zeros!
export update_time_average



"""
    create_ũ_vector(zfv1::AbstractVector)

It allocates the vector to keep in memory the velocity field up to previous 4 time steps
"""
function create_ũ_vector(zfv1::AbstractVector)
    return [deepcopy(zfv1), deepcopy(zfv1), deepcopy(zfv1), deepcopy(zfv1)]
end


"""
    update_ũ_vector!(ũ_vec::Vector, uh_new::AbstractVector)

It updates the vector which stores the values of velocity at previous time steps.
"""
function  update_ũ_vector!(ũ_vec::Vector, uh_new::AbstractVector)
  circshift!(ũ_vec,-1)
  ũ_vec[1] = deepcopy(uh_new)
end

"""
  update_ũ(ũ_vec::Vector)

It uses the Taylor expansion proposed by [Banyai2016](@cite)
"""
function update_ũ(ũ_vec::Vector)
  coeff = [2.1875, -2.1875, 1.3125, -0.3125]
  updt_ũ = ũ_vec[1]*coeff[1] + ũ_vec[2] *coeff[2] + ũ_vec[3] *coeff[3] + ũ_vec[4]*coeff[4]
    return updt_ũ
end



#Extensions for scripts semplification
function GridapDistributed.change_ghost(a::PVector,M::PSparseMatrix)
  col_part = M.col_partition
  GridapDistributed.change_ghost(a,col_part)
end
  
function GridapDistributed.change_ghost(a::PVector,b::AbstractArray)
  GridapDistributed.change_ghost(a,PRange(b))
end

function pazeros(M::PSparseMatrix)
  PartitionedArrays.pzeros(M.col_partition)
end

function  pazeros(a::PVector)
  PartitionedArrays.pzeros(a.index_partition)
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


function update_time_average(field_tn, field_avg, U, tn::Float64, time_step_idx::Int64, time_step::Vector{Float64}, timep::TimeParameters)

  if timep.time_average

    tw0 = timep.time_window[1]
    twf = timep.time_window[2]
    if tw0<=tn<=twf
    
    first_time_step = findfirst(x->x==tw0,time_step)
    last_time_step = findfirst(x->x==twf,time_step)
    n_time_step = findfirst(x->x==tn,time_step)

    n_window_step = n_time_step - first_time_step +1 

    Ntimes = last_time_step-first_time_step+1
    @assert n_window_step<= Ntimes

    @info "updating time average"
  
    field_avg = update_avg_dofs(field_tn, field_avg, U, n_window_step)
  
  end

end

return field_avg

end


function update_avg_dofs(field_tn, field_avg, U, n_step::Int64)
  f_tn_free_dofs = get_free_dof_values(field_tn)

  favg_free_dofs = get_free_dof_values(field_avg)
  favg_free_dofs_update = get_free_dof_values(field_tn)

  if n_step>1
    favg_free_dofs_update = favg_free_dofs*(n_step-1)/n_step + f_tn_free_dofs/n_step
  end

  field_avg = FEFunction(U, favg_free_dofs_update)

  return field_avg
end



end #end module