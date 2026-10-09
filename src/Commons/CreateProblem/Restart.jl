#For periodic case look for ghost nodes, that's the isempty function; The y direction is not periodic for the channel; look for genralization?

"""
    create_search_tree(restart_df::DataFrame)

Create a BruteTree from the NearestNeighbors.jl package taking the coordinates points as input
"""
function create_search_tree(restart_df::DataFrame)
  data = []
  s = 2 #it stores the info if the dataframe is 3D or 2D
  if maximum(abs.(restart_df.uh_2)) > 1e-10
    @info "3D Intialization - it can take a while"
    data = vcat(restart_df.Points_0',restart_df.Points_1',restart_df.Points_2')
    s = 3
  else
    @info "2D Intialization - third velocity component = 0"
    data = vcat(restart_df.Points_0',restart_df.Points_1')
    s = 2
  end

  brutetree = BruteTree(data; leafsize = 12)

  return brutetree, s
end



"""
    restart_uh_field(D::Int64,tree,restart_df::DataFrame)

It provides a suitable function which gives for each point the specified velocity in `restart_file`. It is used as initial condition for restarting a simulation at a specific time step.
"""
function restart_uh_field(D::Int64,s::Int64, tree,restart_df::DataFrame)
  
  if D == 2
    s = 2
  end

  @assert s <= D "Restart DataFrame is $(s) dimensional, the dimension of the problem is $(D). Set uh_2 column == 0 to solve the problem"

  @info "Restarting uh0 ..."
  function u0(x)
    p = [x...][1:s]
    idx, _ = nn(tree, p)
    if D == 2
      return VectorValue(restart_df.uh_0[idx], restart_df.uh_1[idx])
    elseif D==3 && s == 2
        return VectorValue(restart_df.uh_0[idx],  restart_df.uh_1[idx], 0.0)
    elseif D == 3 && s == 3
        return VectorValue(restart_df.uh_0[idx],  restart_df.uh_1[idx], restart_df.uh_2[idx])
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
