#For periodic case look for ghost nodes, that's the isempty function; The y direction is not periodic for the channel; look for genralization?

"""
    find_idx(p::VectorValue{2, Float64}, params)

For a given point `p` it search on the `restart_file` provided by the user the closer node and provides the index of such node.
It takes into account also the channel periodic case, changing the coordiantes of the periodic directions. 
"""

function find_idx(p::VectorValue{2, Float64}, params)
  norm_v = map((x,y)->  norm([p[1]-x, p[2]-y]),  params[:restart_df].Points_0, params[:restart_df].Points_1)
  _,idx = findmin(norm_v)


 
  return idx
end

"""
    find_idx(p::VectorValue{3, Float64}, params)

For a given point `p` it search on the `restart_file` provided by the user the closer node and provides the index of such node.
It takes into account also the channel periodic case, changing the coordiantes of the periodic directions. 
"""
function find_idx(p::VectorValue{3, Float64}, params)
  norm_v = map((x,y,z)->  norm([p[1]-x, p[2]-y,p[3]-z]),  params[:restart_df].Points_0, params[:restart_df].Points_1, params[:restart_df].Points_2)
  _,idx = findmin(norm_v)
  return idx

end



function uh_r(p::VectorValue{2, Float64}, params::Dict{Symbol, Any}, idx::Int)
@unpack restart_df, initial_rescale_factor = params
  VectorValue(restart_df.uh_0[idx][1], restart_df.uh_1[idx][1])*initial_rescale_factor
end


function uh_r(p::VectorValue{3, Float64}, params::Dict{Symbol, Any}, idx::Int)
  @unpack restart_df, initial_rescale_factor = params
  VectorValue(restart_df.uh_0[idx][1], restart_df.uh_1[idx][1],restart_df.uh_2[idx][1])*initial_rescale_factor
end



function ph_r(params::Dict{Symbol, Any}, idx::Int64)
  ph = params[:restart_df].ph[idx][1]
    return  ph
end



"""
    uh_restart(p, params::Dict{Symbol, Any})

For a given point `p` it calles `find_idx` which provide the line of the `csv` file corresponding to that point. Then, it calles `uh` which provide the `VectorValue` of the velocity at that point.
"""
function uh_restart(p, params::Dict{Symbol, Any})
  
  idx = find_idx(p, params)
  return uh_r(p, params, idx)

end


"""
    ph_restart(p, params::Dict{Symbol, Any})

For a given point `p` it calles `find_idx` which provide the line of the `csv` file corresponding to that point. Then, it calles `ph` which provide the scalar pressure value at that point.
"""
function ph_restart(p, params::Dict{Symbol, Any})
  idx = find_idx(p, params)
  ph = ph_r(params,idx)
  return ph
end

"""
    restart_uh_field(params::Dict{Symbol, Any})

It provides a suitable function which gives for each point the specified velocity in `restart_file`. It is used as initial condition for restarting a simulation at a specific time step.
"""
function restart_uh_field(params::Dict{Symbol, Any})
  println("Restarting uh0 ...")
  # u0(t::Real) = x -> u0(x, t::Real)
  u0(x) = uh_restart(x, params)

  return u0

end

"""
    restart_ph_field(params::Dict{Symbol, Any})

It provides a suitable function which gives for each point the specified pressure in `restart_file`. It is used as initial condition for restarting a simulation at a specific time step.
"""
function restart_ph_field(params::Dict{Symbol, Any})
  println("Restarting ph0 ...")
  
  init_pres = DataFrames.haskey(params[:restart_df],:ph)

  # p0(x, t::Real) = (init_pres) ?   ph_restart(x, params) : 0.0
  # p0(t::Real) = x -> p0(x, t::Real)
  p0(x) = (init_pres) ?   ph_restart(x, params) : 0.0

  return p0

end
