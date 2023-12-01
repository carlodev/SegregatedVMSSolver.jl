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
export get_local_unique_idx
export export_nodes_glob
export export_n_Γ

export create_export_tags!



function unwrap_vector(a)
    V = Float64[]
for v in a
    append!(V,v...)
end    
return V
end

function wrap_vector(vv,p::Int64)
    L = length(vv)
    l = Int(L/p)
    V = [Vector{Float64}(undef,p) for _ in 1:l]
for i in 1:l
    for j in 1:p
        idx = (i-1)*p + j
        V[i][j] = vv[idx]
    end
end    

return V
end




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
    D = get_dimension(vv)
    x = getindex.(vv,1)
    y = getindex.(vv,2)
    z = (D==2) ? zeros(length(x)) : getindex.(vv,3)
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
    map(gfield, parts) do g, part
        if part == 1 #On Main procs
            values_uniques = g.data[global_unique_idx]
            export_time_step(timestep, values_uniques, fieldname)
        end
    end
end


function get_visgrid(ttrian)
    f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)
    ref_grids = map(f, Gridap.Geometry.get_reffes(ttrian))
    visgrid = Gridap.Visualization.VisualizationGrid(ttrian, ref_grids)
    return visgrid
end

# """
#     get_local_unique_idx(parts, trian)

# For each part (aka processors), only non-duplicated nodes indexes are extracted.
# """
function get_local_unique_idx(params)

    @unpack parts, export_tags = params
    local_unique_idx = nothing

    if !isnothing(export_tags)
        @unpack Γ_ = export_tags
        local_unique_idx = []
        for (Γ) in zip(Γ_)

            #export nodes
            local_unique_idx_ = map(Γ.trians) do ttrian
                visgrid = get_visgrid(ttrian)
                visgrid_ = conv_VectorValue.(visgrid.sub_grid.node_coordinates)
                unique_idx = unique(i -> visgrid_[i], eachindex(visgrid_)) #Indexes of unique nodes on each part
                return unique_idx
            end
            push!(local_unique_idx,local_unique_idx_)
        end #end for

        merge!(export_tags, Dict(:local_unique_idx_=>local_unique_idx))
    end #end if 
    return local_unique_idx
end




# """
#     export_nodes_glob(parts, trian)

# It gathers on MAIN procs (==1) the non-duplicates nodes for each procs. It computes the non-duplicate global indexes.
# It also extracts non-duplicated nodes.
# """
function export_nodes_glob(parts, trian, tagname,D)

    #export nodes
    local_unique_nodes = map(trian.trians) do ttrian
        visgrid = get_visgrid(ttrian)
        visgrid_ = conv_VectorValue.(visgrid.sub_grid.node_coordinates)
        nodes_tri = unique(visgrid_)
        return unwrap_vector(nodes_tri)
    end

    
    glun = gather(local_unique_nodes)

    global_unique_idx = 0.0
    map(glun,parts) do g,part
        if part ==1 #On Main procs
            gg =  wrap_vector(g.data,D)
            nodes_uniques = unique(gg)

            export_time_step(0.0, nodes_uniques, "$(tagname)_nodes")
            global_unique_idx = unique(i -> g.data[i], eachindex(gg)) 
        end
    end
    return global_unique_idx
end


function export_nodes_glob(params::Dict{Symbol,Any})
    @unpack parts, export_tags, model,D = params
  
    global_unique_idx = nothing

    if !isnothing(export_tags)
        global_unique_idx = []
        @unpack name_tags, Γ_ = export_tags
        for (name_tag, Γ) in zip(name_tags, Γ_)
            push!(global_unique_idx,export_nodes_glob(parts, Γ,name_tag,D))
        end
        merge!(export_tags, Dict(:global_unique_idx_=>global_unique_idx))
    end
end


# """
#     export_n_Γ(params::Dict{Symbol,Any}, local_unique_idx, global_unique_idx)

# Export tangent and normal vectors at the corresponding points.
# """
function export_n_Γ(params::Dict{Symbol,Any})
    #get unique values in each processor

    @unpack export_tags,parts = params


    if !isnothing(export_tags)
        @unpack name_tags, Γ_, n_Γ_, local_unique_idx_, global_unique_idx_ = export_tags

        for (name_tag, Γ, n_Γ,local_unique_idx, global_unique_idx) in zip(name_tags, Γ_, n_Γ_,local_unique_idx_, global_unique_idx_)
           
            n_Γ = -n_Γ #pointing from the body to the outside
            t_Γ = rotation ∘ n_Γ #extract tangent

            cellfields = Dict("$(name_tag)_n_Γ" => n_Γ, "$(name_tag)_t_Γ" => t_Γ)

            fdat = GridapDistributed._prepare_fdata(Γ.trians, cellfields)
            for field in keys(cellfields)
                fieldh = map(Γ.trians, fdat, local_unique_idx) do ttrian, cf, unique_idx
                    visgrid = get_visgrid(ttrian)
                    pdata = Gridap.Visualization._prepare_pdata(ttrian, cf, visgrid.cell_to_refpoints)
                    field_h = 0.0

                    field_h = pdata[field][unique_idx]
                    return field_h
                end
                extract_global_unique(fieldh, parts, global_unique_idx, 0.0, field)
            end

        end #for

    end #end if
end

"""
    export_fields(params::Dict{Symbol,Any}, local_unique_idx, global_unique_idx, tt::Float64, uh0, ph0)

Export pressure and friction (not multiplied by viscosity) - airfoil simulations oriented
"""
function export_fields(params::Dict{Symbol,Any},  tt::Float64, uh0, ph0)
    #get unique values in each processor

    @unpack parts, export_tags = params

    if !isnothing(export_tags)
        @unpack name_tags, Γ_, n_Γ_, local_unique_idx_, global_unique_idx_ = export_tags
        @unpack fieldexport = params
        for (name_tag, Γ, n_Γ,local_unique_idx, global_unique_idx,fieldexp) in zip(name_tags, Γ_, n_Γ_,local_unique_idx_, global_unique_idx_,fieldexport)
            n_Γ = -n_Γ #pointing from the body to the outside
            t_Γ = rotation ∘ n_Γ #extract tangent

            friction = (transpose(∇(uh0)) ⋅ n_Γ) ⋅ t_Γ

            cellfields = create_cellfield(name_tag,uh0,ph0,friction,fieldexp)
            fdat = GridapDistributed._prepare_fdata(Γ.trians, cellfields)

            for field in keys(cellfields)
                fieldh = map(Γ.trians, fdat, local_unique_idx) do ttrian, cf, unique_idx
                    visgrid = get_visgrid(ttrian)
                    pdata = Gridap.Visualization._prepare_pdata(ttrian, cf, visgrid.cell_to_refpoints)
                    field_h = pdata[field][unique_idx]
                    return field_h
                end #end map
                extract_global_unique(fieldh, parts, global_unique_idx, tt, field)
            end #end for field in keys(cellfields)

        end #for () in zip ()
    end #end if
end


function create_cellfield(name_tag::String,uh,ph,friction,fieldexp)
    cellfields = Dict()
    for fe in fieldexp
        if fe == "uh"
            cf = Dict("$(name_tag)_uh" => uh)
        elseif fe=="ph"
            cf = Dict("$(name_tag)_ph" => ph)
        elseif fe=="friction"
            cf = Dict("$(name_tag)_friction" => friction)
        else
            @error "export field $fe not recognized as valid, use uh,ph,friction"
        end
        merge!(cellfields,cf)
    end
    return cellfields
end


# """
#     rotation(n::VectorValue{2,Float64})

# It rotates by π/2 the n vector
# """
function rotation(n::VectorValue{2,Float64})
    n1, n2 = [n...]
    VectorValue(-n2, n1)
end

function rotation(n::VectorValue{3,Float64})
    n1, n2, n3 = [n...]
    VectorValue(-n2, n1, n3)
end



"""
    create_export_tags!(params::Dict{Symbol,Any})

For each ´name_tags´ it creates the ´export_tags´ dictionary
"""
function create_export_tags!(params::Dict{Symbol,Any})

    @unpack model,name_tags = params
    export_tags = nothing # If no name tags are specified to export
    if !isnothing(name_tags)
        Γ = []
        n_Γ = []
        for nt in name_tags
            Γ_tmp = BoundaryTriangulation(model; tags=nt)
            n_Γ_tmp = get_normal_vector(Γ_tmp)
            push!(Γ,Γ_tmp)
            push!(n_Γ,n_Γ_tmp)
        end
        export_tags = Dict(:name_tags => name_tags, :Γ_ => Γ, :n_Γ_ => n_Γ)
    end
    merge!(params, Dict(:export_tags=>export_tags))
end

"""
    print_on_request(log_dir::String)

If in the directory `log_dir` exists and there is the `PrintSim.txt` file return true. If true is printing simulations results in Paraview format (.vtu and .pvtu). 
"""
function print_on_request(log_dir::String)
    flag = false
    current_dir = readdir()
    log_idx = findfirst(x -> x == log_dir, current_dir)

    if !isnothing(log_idx)
        log_dir = readdir(current_dir[log_idx])
        print_idx = findfirst(x -> x == "PrintSim.txt", log_dir)
        if !isnothing(print_idx)
            flag = true
        end
    end

    return flag
end

end



