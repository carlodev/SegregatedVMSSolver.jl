
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


function Base.println(d::Dict)
  for k in keys(d)
    if k !== :restart_df
    kval = d[k]
    println("$k = $kval")
    end
  end
end


function DataFrames.haskey(df::DataFrame,s::Symbol)
  
  if sum(DataFrames.propertynames(df) .== s)>0
        return true
  else
        return false
  end

end

function update_time_average!(uh_tn,ph_tn, uh_avg, ph_avg, tn, params)
  @unpack time_step, time_window= params
  if !isnothing(time_window)
    tw0 = time_window[1]
    twf = time_window[2]
    if tw0<=tn<=twf
    
    tsc=collect(time_step)
    first_time_step = findfirst(x->x==tw0,tsc)
    last_time_step = findfirst(x->x==twf,tsc)
    
    Ntimes = last_time_step-first_time_step+1
    println("updating time average")

    update_avg_dofs!(deepcopy(uh_tn), uh_avg, Ntimes)
    update_avg_dofs!(deepcopy(ph_tn), ph_avg, Ntimes)
  end

end
end




function update_avg_dofs!(f_tn, f_avg, N)
  println("updates fields")
  update_avg_dofs!(f_tn.fields,  f_avg.fields, N)
end

function update_avg_dofs!(fields_tn::MPIArray, fields_avg::MPIArray, N)
  println("updates fields")
  fields_avg.item.free_values .+= fields_tn.item.free_values ./N
end


function update_avg_dofs!(fields_tn::DebugArray, fields_avg::DebugArray, N)
  println("updates fields")
  for (ff_tn, ff_avg) in zip(fields_tn,fields_avg)
    ff_avg.free_values .+= ff_tn.free_values ./N
  end

end