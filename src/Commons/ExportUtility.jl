module ExportUtility

using Gridap
using GridapDistributed
using CSV, DelimitedFiles
using DataFrames
using Parameters
using PartitionedArrays

export conv_to_df
export writesolution
export initialize_export_nodes
export export_fields

using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.Interfaces

function unwrap_vector(a)
    V = Float64[]
    for v in a
        append!(V, v...)
    end
    return V
end

function wrap_vector(vv, p::Int64)
    L = length(vv)
    l = Int(L / p)
    V = [Vector{Float64}(undef, p) for _ in 1:l]
    for i in 1:l
        for j in 1:p
            idx = (i - 1) * p + j
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
    x = getindex.(vv, 1)
    y = getindex.(vv, 2)
    z = (D == 2) ? zeros(length(x)) : getindex.(vv, 3)
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
#     export_time_step(t::Float64, vv::Vector, fname::String)

# Save .csv file global nodes
# """
function export_time_step(t::Float64, vv::Vector, fname::String)
    df = conv_to_df(vv)
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
    @unpack Γ_ = export_tags

    local_unique_idx = []
    for Γ in Γ_

        #export nodes
        local_unique_idx_ = map(Γ.trians) do ttrian
            visgrid = get_visgrid(ttrian)
            visgrid_ = conv_VectorValue.(visgrid.sub_grid.node_coordinates)
            unique_idx = unique(i -> visgrid_[i], eachindex(visgrid_)) #Indexes of unique nodes on each part
            return unique_idx
        end
        push!(local_unique_idx, local_unique_idx_)
    end #end for

    merge!(export_tags, Dict(:local_unique_idx_ => local_unique_idx))

    return local_unique_idx
end




# """
#     export_nodes_glob(parts, trian)

# It gathers on MAIN procs (==1) the non-duplicates nodes for each procs. It computes the non-duplicate global indexes.
# It also extracts non-duplicated nodes.
# """
function export_nodes_glob(parts, trian, tagname::String, D::Int64)

    #export nodes
    local_unique_nodes = map(trian.trians) do ttrian
        visgrid = get_visgrid(ttrian)
        visgrid_ = conv_VectorValue.(visgrid.sub_grid.node_coordinates)
        nodes_tri = unique(visgrid_)
        return unwrap_vector(nodes_tri)
    end


    glun = gather(local_unique_nodes)

    global_unique_idx = 0.0
    map(glun, parts) do g, part
        if part == 1 #On Main procs
            gg = wrap_vector(g.data, D)
            nodes_uniques = unique(gg)

            export_time_step(0.0, nodes_uniques, "$(tagname)_nodes")
            global_unique_idx = unique(i -> gg[i], eachindex(gg))
        end
    end
    return global_unique_idx
end


function export_nodes_glob(params::Dict{Symbol,Any}, D::Int64)
    @unpack parts, export_tags, model = params
    @unpack name_tags, Γ_ = export_tags

    global_unique_idx = []
    for (name_tag, Γ) in zip(name_tags, Γ_)
        push!(global_unique_idx, export_nodes_glob(parts, Γ, name_tag, D))
    end
    merge!(export_tags, Dict(:global_unique_idx_ => global_unique_idx))
end


# """
#     export_n_Γ(params::Dict{Symbol,Any}, local_unique_idx, global_unique_idx)

# Export tangent and normal vectors at the corresponding points.
# """
function export_n_Γ(params::Dict{Symbol,Any})
    #get unique values in each processor

    @unpack export_tags, parts = params


    @unpack name_tags, Γ_, n_Γ_, local_unique_idx_, global_unique_idx_ = export_tags

    for (name_tag, Γ, n_Γ, local_unique_idx, global_unique_idx) in zip(name_tags, Γ_, n_Γ_, local_unique_idx_, global_unique_idx_)

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

end


function export_fields(params, simcase, tt, uh0, ph0)
    @sunpack export_field, fieldexport = simcase.simulationp.exportp
    if export_field
        export_fields(params::Dict{Symbol,Any}, fieldexport, tt::Float64, uh0, ph0)
    end
end


"""
    export_fields(params::Dict{Symbol,Any}, fieldexport::Vector, tt::Float64, uh0, ph0)
Export pressure and friction (not multiplied by viscosity) - airfoil simulations oriented
"""
function export_fields(params::Dict{Symbol,Any}, fieldexport::Vector, tt::Float64, uh0, ph0)
    #get unique values in each processor

    @unpack parts, export_tags = params
    @unpack name_tags, Γ_, n_Γ_, local_unique_idx_, global_unique_idx_ = export_tags

    for (name_tag, Γ, n_Γ, local_unique_idx, global_unique_idx, fieldexp) in zip(name_tags, Γ_, n_Γ_, local_unique_idx_, global_unique_idx_, fieldexport)
        n_Γ = -n_Γ #pointing from the body to the outside
        t_Γ = rotation ∘ n_Γ #extract tangent

        friction = (transpose(∇(uh0)) ⋅ n_Γ) ⋅ t_Γ

        cellfields = create_cellfield(name_tag, uh0, ph0, friction, fieldexp)
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
end


function create_cellfield(name_tag::String, uh, ph, friction, fieldexp)
    cellfields = Dict()
    for fe in fieldexp
        if fe == "uh"
            cf = Dict("$(name_tag)_uh" => uh)
        elseif fe == "ph"
            cf = Dict("$(name_tag)_ph" => ph)
        elseif fe == "friction"
            cf = Dict("$(name_tag)_friction" => friction)
        else
            @error "export field $fe not recognized as valid, use uh,ph,friction"
        end
        merge!(cellfields, cf)
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

function initialize_export_nodes(params::Dict{Symbol,Any}, simcase::SimulationCase)
    @sunpack export_field, name_tags = simcase.simulationp.exportp
    @sunpack D = simcase
    if export_field
        mkpath("Results")
        @info "Folder Results created"
        
        create_export_tags!(params, name_tags)
        get_local_unique_idx(params)
        export_nodes_glob(params, D)
        export_n_Γ(params)
        @info "Nodes Coordinates Exported"
    end

end

"""
    create_export_tags!(params::Dict{Symbol,Any})

For each ´name_tags´ it creates the ´export_tags´ dictionary
"""
function create_export_tags!(params::Dict{Symbol,Any}, name_tags::Vector)
    @unpack model = params

    export_tags = Dict(:name_tags => name_tags)
    Γ = []
    n_Γ = []
    for nt in name_tags
        Γ_tmp = BoundaryTriangulation(model; tags=nt)
        n_Γ_tmp = get_normal_vector(Γ_tmp)
        push!(Γ, Γ_tmp)
        push!(n_Γ, n_Γ_tmp)
    end
    export_tags = Dict(:name_tags => name_tags, :Γ_ => Γ, :n_Γ_ => n_Γ)

    merge!(params, Dict(:export_tags => export_tags))
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


function writesolution(params::Dict{Symbol,Any}, simcase::SimulationCase, ntime::Int64, tn::Float64, fields::Tuple, avg_fields::Tuple)
    benchmark = simcase.simulationp.exportp.benchmark
    log_dir = simcase.simulationp.exportp.log_dir
    nsub = simcase.sprob.method.order

    if !benchmark
        if (mod(ntime, 100) == 0 || print_on_request(log_dir))
            save_sim_dir = simcase.simulationp.exportp.save_sim_dir
            case = typeof(simcase)
            save_path = joinpath(save_sim_dir, "$(case)_$(tn)_.vtu")
        
            writesolution(simcase, params, nsub, save_path, tn, fields,avg_fields)
        end
        compute_enstrophy_ek(simcase, params,fields,tn)

    end
    compute_error(simcase, params,   tn, fields)

end


function writesolution(simcase, params::Dict{Symbol,Any}, nsub::Int64, save_path::String, tn::Float64, fields::Tuple, avg_fields::Tuple)
    uh, ph = fields
    uh_avg, ph_avg =avg_fields
    fields_names = simcase.simulationp.exportp.vtu_export
    @unpack Ω = params

    ALLOWED_EXPORTS = ["uh", "ph", "uh_avg", "ph_avg"]
    selected_exports = intersect(fields_names,ALLOWED_EXPORTS)
    cell_fields = Dict{String,Any}()
    if "uh" in selected_exports
        merge!(cell_fields,Dict("uh"=>uh))
    end
    if "ph" in selected_exports
        merge!(cell_fields,Dict("ph"=>ph))
    end
    if "uh_avg" in selected_exports
        merge!(cell_fields,Dict("uh_avg"=>uh_avg))
    end
    if "ph_avg" in selected_exports
        merge!(cell_fields,Dict("ph_avg"=>ph_avg))
    end
   
    writesolution(simcase, tn, fields_names,cell_fields)

    writevtk(Ω, save_path, nsubcells=nsub, cellfields=cell_fields)
end

function writesolution(simcase, tn::Float64, fields_names, cell_fields::Dict{String,Any})
end

function writesolution(simcase::TaylorGreen{Periodic}, tn::Float64, fields_names, cell_fields::Dict{String,Any})
    @sunpack D = simcase

    if D == 2
        ALLOWED_EXPORTS = ["uh_analytic", "ph_analytic"]
            uh_analytic = simcase.bc_type.a_solution[:velocity](tn) # Compute analytic velocity at tn
            ph_analytic = simcase.bc_type.a_solution[:pressure](tn) # Compute analytic pressure at tn
            selected_exports = intersect(fields_names,ALLOWED_EXPORTS)


            if "uh_analytic" in selected_exports
                merge!(cell_fields,Dict("uh_analytic"=>uh_analytic))
        
            end
            if "ph_analytic" in selected_exports
                merge!(cell_fields,Dict("ph_analytic"=>ph_analytic))
            end


    end

end


"""
    compute_error(params::Dict{Symbol,Any}, simcase::TaylorGreen{Periodic},  tn::Float64, fields::Tuple)

Compute L2 error norm for velocity and pressure for TGV2D case. https://doi.org/10.1016/j.enganabound.2020.12.018
"""
function compute_error(simcase::TaylorGreen{Periodic}, params::Dict{Symbol,Any},   tn::Float64, fields::Tuple)
    uh, ph = fields

    uh_analytic = simcase.bc_type.a_solution[:velocity](tn) # Compute analytic velocity at tn
    ph_analytic = simcase.bc_type.a_solution[:pressure](tn) # Compute analytic pressure at tn

    @unpack dΩ,parts = params
    eu = uh_analytic - uh    #error velocity 
    ep = ph_analytic - ph #error pressure

    #L2 norm error velocity and pressure
    l2eu = sqrt(sum(∫(eu ⋅ eu) * dΩ)) ./sqrt(sum(∫(uh ⋅ uh) * dΩ))
    l2ep = sqrt(sum(∫(ep * ep) * dΩ)) ./ sqrt(sum(∫(ph ⋅ ph) * dΩ))

    ALLOWED_EXPORTS = ["VelocityError", "PressureError"]
    selected_exports = intersect(ALLOWED_EXPORTS, simcase.simulationp.exportp.extra_export)

    # Construct data array dynamically
    data = [tn]
    headers = ["time"]
    if "VelocityError" in selected_exports
        push!(data, l2eu)
        push!(headers, "VelocityError")
    end
    if "PressureError" in selected_exports
        push!(data, l2ep)
        push!(headers, "PressureError")

    end

    write_to_csv("TGV_ERRRORS.csv", data, headers, parts)
    
end

function compute_error(simcase, params,   tn, fields)
    
end

"""
    compute_enstrophy_ek(simcase, params::Dict{Symbol,Any}, fields::Tuple)

Compute enstrophy and Kinetic Energy - not normalized with volume.
"""
function compute_enstrophy_ek(simcase, params::Dict{Symbol,Any}, fields::Tuple, tn::Float64)
    uh, _ = fields
    @unpack dΩ,parts = params

    wh = ∇×uh

    ### Compute Kinetic Energy
    Ek = 0.5 .* sum(∫(uh ⋅ uh) * dΩ) #not normalized with the volume

    ### Compute Dissipation or Enstrophy https://en.wikipedia.org/wiki/Enstrophy
    Enstrophy = sum(∫(wh ⋅ wh) * dΩ) #not normalized with the volume

    

    ALLOWED_EXPORTS = ["Enstrophy", "KineticEnergy"]
    selected_exports = intersect(ALLOWED_EXPORTS, simcase.simulationp.exportp.extra_export)


    # Construct data array dynamically
    data = [tn]
    headers = ["time"]
    if "Enstrophy" in selected_exports
        push!(data, Enstrophy)
        push!(headers, "Enstrophy")

    end
    if "KineticEnergy" in selected_exports
        push!(data, Ek)
        push!(headers, "KineticEnergy")

    end

    write_to_csv("Enstrophy_EK.csv", data, headers, parts)

end



"""
    write_to_csv(file_path::String, data::Vector{Vector{Any}}, parts)

Write Results on CSV file
"""
function write_to_csv(file_path::String, data::Vector{Float64}, headers::Vector{String}, parts)
    if length(data) > 1 #if length(data) ==1 only would be time exported
    map( parts) do part
        if part == 1
            
            file_exists = isfile(file_path)

            open(file_path, "a") do file
                # Write headers if the file does not exist
                if !file_exists
                    writedlm(file, [headers], ',')
                end
            end

            open(file_path, "a") do file
                writedlm(file, [data], ',')
            end
        end
    end
end

end




end # end module