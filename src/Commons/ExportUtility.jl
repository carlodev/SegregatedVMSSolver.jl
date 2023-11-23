module ExportUtility

using Gridap
using GridapDistributed
using CSV
using DataFrames
using Parameters
using PartitionedArrays

export conv_to_df
export print_on_request
export export_fields
export forces_domain!
export get_local_unique_idx
export export_nodes_glob
export export_n_Γ

# """
#     conv_VectorValue(v::VectorValue)

# Convert VectorValue (Gridap Type) to a Vector
# """
function conv_VectorValue(v::VectorValue)
    [v...]
end

# """
#     get_dimension(vv::Vector)

# It provides the dimension
# """
function get_dimension(vv::Vector)
    D = 0
    if length(vv) > 0
        D = length(vv[1])
    end
    return D

end

"""
    conv_to_df(vv::Vector)

Convert a Vector{Vector} to a DataFrame. It is used for export nodes, normals.
"""
function conv_to_df(vv::Vector)
    d = get_dimension(vv)
    n = length(vv)
    x = zeros(n)
    y = zeros(n)
    z = zeros(n)
    for (i, val) in enumerate(vv)
        x[i] = val[1]
        y[i] = val[2]
        if d == 3
            z[i] = val[3]
        end

    end
    df = DataFrame(x=x, y=y, z=z)
    return df
end


"""
    conv_to_df(vv::Vector)

Convert a Vector{Float64} to a DataFrame. It is used for export pressure field.
"""
function conv_to_df(vv::Vector{Float64})
    df = DataFrame(p=vv)
    return df
end

# """
#     export_time_step(t::Float64, vv::Vector, fname::String, part::Int64)

# Save .csv file, one for each processor
# """
function export_time_step(t::Float64, vv::Vector, fname::String, part::Int64)
    df = conv_to_df(vv)

    dir = joinpath("Results", "$(fname)_$t")
    mkpath(dir)
    filename = joinpath(dir, "$(fname)_$(part)_$t.csv")
    CSV.write(filename, df)
end

# """
#     export_time_step(t::Float64, vv::Vector, fname::String)

# Save .csv file global nodes
# """
function export_time_step(t::Float64, vv::Vector, fname::String)
    df = conv_to_df(vv)
    mkpath("Results")
    filename = joinpath("Results", "$(fname)_$t.csv")
    CSV.write(filename, df)
end


function extract_global_unique(dfield, parts, global_unique_idx, timestep::Float64, fieldname::String)
    gfield = gather(dfield)
    map(gfield,parts) do g, part
        if part ==1 #On Main procs
            values_uniques = g.data[global_unique_idx]
            export_time_step(timestep, values_uniques, fieldname)          
        end
    end
end


# """
#     get_local_unique_idx(parts, trian)

# For each part (aka processors), only non-duplicated nodes indexes are extracted.
# """
function get_local_unique_idx(parts, trian)
    f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)

    #export nodes
    local_unique_idx = map(parts, trian.trians) do part, ttrian
        ref_grids = map(f, Gridap.Geometry.get_reffes(ttrian))
        visgrid = Gridap.Visualization.VisualizationGrid(ttrian, ref_grids)
        visgrid_ = conv_VectorValue.(visgrid.sub_grid.node_coordinates)
        nodes_tri = unique(visgrid_) #Coordinate of unique nodes
 
        unique_idx = unique(i -> visgrid_[i], eachindex(visgrid_)) #Indexes of unique nodes on each part
        # export_time_step(0.0, nodes_tri, "nodes", part)

        return unique_idx
    end

    return local_unique_idx
end




# """
#     export_nodes_glob(parts, trian)

# It gathers on MAIN procs (==1) the non-duplicates nodes for each procs. It computes the non-duplicate global indexes.
# It also extracts non-duplicated nodes.
# """
function export_nodes_glob(parts, trian)
    f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)

    #export nodes
    local_unique_nodes = map(parts, trian.trians) do part, ttrian
        ref_grids = map(f, Gridap.Geometry.get_reffes(ttrian))
        visgrid = Gridap.Visualization.VisualizationGrid(ttrian, ref_grids)
        visgrid_ = conv_VectorValue.(visgrid.sub_grid.node_coordinates)
        nodes_tri = unique(visgrid_)
        return nodes_tri
    end

    glun = gather(local_unique_nodes)

    global_unique_idx = 0.0
    map(glun,parts) do g,part
        if part ==1 #On Main procs

            nodes_uniques = unique(g.data)
            export_time_step(0.0, nodes_uniques, "nodes")

            global_unique_idx = unique(i -> g.data[i], eachindex(g.data)) 
        end
    end
    return global_unique_idx
end




# """
#     export_n_Γ(params::Dict{Symbol,Any}, local_unique_idx, global_unique_idx)

# Export tangent and normal vectors at the corresponding points.
# """
function export_n_Γ(params::Dict{Symbol,Any}, local_unique_idx, global_unique_idx)
    #get unique values in each processor

    @unpack Γ,n_Γ, parts,force_tags = params

    if force_tags !== nothing
        f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)

        n_Γ = - n_Γ #pointing from the body to the outside
        t_Γ = rotation∘n_Γ #extract tangent

        cellfields = Dict( "n_Γ" => n_Γ, "t_Γ"=>t_Γ)

        fdat = GridapDistributed._prepare_fdata(Γ.trians, cellfields)
        for field in keys(cellfields)
        fieldh = map(parts, Γ.trians, fdat, local_unique_idx) do part, ttrian, cf, unique_idx
            ref_grids = map(f, Gridap.Geometry.get_reffes(ttrian))
            visgrid = Gridap.Visualization.VisualizationGrid(ttrian, ref_grids)
            pdata = Gridap.Visualization._prepare_pdata(ttrian, cf, visgrid.cell_to_refpoints)
            field_h = 0.0
            
            field_h = pdata[field][unique_idx]
            # export_time_step(0.0, field_h, field, part)
            return field_h
            end
         extract_global_unique(fieldh, parts, global_unique_idx, 0.0, field)
    end #for
    
    end #end if
end

"""
    export_fields(params::Dict{Symbol,Any}, local_unique_idx, global_unique_idx, tt::Float64, uh0, ph0)

Export pressure and friction (not multiplied by viscosity) - airfoil simulations oriented
"""
function export_fields(params::Dict{Symbol,Any}, local_unique_idx, global_unique_idx, tt::Float64, uh0, ph0)
    #get unique values in each processor

    @unpack Γ,n_Γ, parts,force_tags = params

    if force_tags !== nothing
        f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)

        n_Γ = - n_Γ #pointing from the body to the outside
        t_Γ = rotation∘n_Γ #extract tangent

        friction = (transpose(∇(uh0))⋅n_Γ) ⋅ t_Γ

        cellfields = Dict("ph" => ph0, "friction" => friction)

        fdat = GridapDistributed._prepare_fdata(Γ.trians, cellfields)

        for field in keys(cellfields)
            fieldh = map(parts, Γ.trians, fdat, local_unique_idx) do part, ttrian, cf, unique_idx
                ref_grids = map(f, Gridap.Geometry.get_reffes(ttrian))
                visgrid = Gridap.Visualization.VisualizationGrid(ttrian, ref_grids)
                pdata = Gridap.Visualization._prepare_pdata(ttrian, cf, visgrid.cell_to_refpoints)
                field_h = pdata[field][unique_idx]
            return field_h 
            end #end map
            extract_global_unique(fieldh, parts, global_unique_idx, tt, field)
        end #end for
    end
end

# """
#     rotation(n::VectorValue{2,Float64})

# It rotates by π/2 the n vector
# """
function rotation(n::VectorValue{2,Float64})
    n1,n2 = [n...]
    VectorValue(-n2,n1)
end

function rotation(n::VectorValue{3,Float64})
    n1,n2,n3 = [n...]
    VectorValue(-n2,n1,n3)
end


# """
#     forces_domain(model, force_tags, degree)

# For a given `force_tags` in the `model` it provides the triangulation, the measure and the normal vector.
# """
function forces_domain!(params)
    @unpack model,force_tags = params
    degree = 4

    Γ = BoundaryTriangulation(model; tags=force_tags) 
    dΓ = Measure(Γ,degree)
    n_Γ = get_normal_vector(Γ)
    force_params = Dict(:Γ => Γ, :dΓ => dΓ, :n_Γ => n_Γ)
    merge!(params,force_params)

end


"""
    print_on_request(log_dir::String)

If in the directory `log_dir` exists and there is the `PrintSim.txt` file return true. If true is printing simulations results in Paraview format (.vtu and .pvtu). 
"""
function print_on_request(log_dir::String)
    flag = false
    current_dir = readdir()
    log_idx = findfirst(x->x==log_dir, current_dir)
    
    if !isnothing(log_idx)
        log_dir = readdir(current_dir[log_idx])
        print_idx = findfirst(x->x =="PrintSim.txt", log_dir)
        if !isnothing(print_idx)
            flag = true
        end
    end
    
    return flag
end

end