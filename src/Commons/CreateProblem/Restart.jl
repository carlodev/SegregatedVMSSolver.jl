#For periodic case look for ghost nodes, that's the isempty function; The y direction is not periodic for the channel; look for genralization?

"""
    create_search_tree(restart_df::DataFrame)

Create a BruteTree from the NearestNeighbors.jl package taking the coordinates points as input
"""
function create_search_tree(restart_df::DataFrame)

  data = vcat(restart_df.Points_0',restart_df.Points_1')
  brutetree = BruteTree(data; leafsize = 12)

  return brutetree
end



"""
    restart_uh_field(D::Int64,tree,restart_df::DataFrame)

It provides a suitable function which gives for each point the specified velocity in `restart_file`. It is used as initial condition for restarting a simulation at a specific time step.
"""
function restart_uh_field(D::Int64,tree,restart_df::DataFrame)

  @info "Restarting uh0 ..."

  function u0(x)
    p = [x...][1:2]
    idx, _ = nn(tree, p)
    if D == 2
      return VectorValue(restart_df.uh_0[idx], restart_df.uh_1[idx])
    elseif D==3
      return VectorValue(restart_df.uh_0[idx],  restart_df.uh_1[idx], 0.0)
    end
  end

  return u0

end

"""
    restart_ph_field(simcase::VelocityBoundary,tree)

It provides a suitable function which gives for each point the specified pressure in `restart_file`. It is used as initial condition for restarting a simulation at a specific time step.
"""
function restart_ph_field(tree,restart_df::DataFrame)

  @info "Restarting ph0 ..."
  
  init_pres = DataFrames.haskey(restart_df,:ph)

  function ph_restart(x)
    p = [x...][1:2]
    idx, _ = nn(tree, p)
    ph = restart_df.ph[idx]
    return ph
  end

  p0(x) = (init_pres) ?   ph_restart(x) : 0.0

  return p0

end

function DataFrames.haskey(df::DataFrame,s::Symbol)
  
  if sum(DataFrames.propertynames(df) .== s)>0
        return true
  else
        return false
  end

end
