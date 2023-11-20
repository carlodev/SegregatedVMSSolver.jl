
function custom_cmp_file(x::String)
    str_len = length(x)

    offset = findlast("/",x)[1]
    number_idx = findfirst(isdigit, x[offset:end]) + offset-1
    str, num = SubString(x, 1, number_idx - 1), SubString(x, number_idx, str_len-4)
    return str, parse(Float64, num)
end

function custom_cmp_dir(x::String)
    str_len = length(x)

    offset = findlast("/",x)[1]
    number_idx = findfirst(isdigit, x[offset:end]) + offset-1
    str, num = SubString(x, 1, number_idx - 1), SubString(x, number_idx, str_len)
    
    return str, parse(Float64, num)
end


function get_file_dir_idx(path0::String)
    s = readdir(path0)
    path = path0 .* s

    file_idx = Int64[]
    dir_idx = Int64[]
    sdir = isdir.(path)

    for (i, val) in enumerate(sdir)
        if val
            push!(dir_idx, i)
        else
            push!(file_idx, i)
        end
    end

    return file_idx, dir_idx, path

end



function get_file_dir_idx_fields(path::Vector{String}, file_idx::Vector, dir_idx::Vector, cellfields::Dict)
    field_dir_idx = Vector[]

    field_file_idx = Vector[]


    for field_name in keys(cellfields)
        v_dir = map(x -> occursin(field_name, x), path[dir_idx])
        v_file = map(x -> occursin(field_name, x), path[file_idx])

        idx_f_dir = findall(x -> x == true, v_dir)
        idx_f_file = findall(x -> x == true, v_file)

        idx_dir = sortperm(path[dir_idx][idx_f_dir], by=custom_cmp_dir)
        idx_file = sortperm(path[file_idx][idx_f_file], by=custom_cmp_file)

        push!(field_dir_idx, idx_f_dir[idx_dir])
        push!(field_file_idx, idx_f_file[idx_file])

    end

    return field_file_idx, field_dir_idx
end

function move_files(path, dir_idx::Vector)
    isairfoil = occursin.("Airfoil", path[dir_idx])
    airfoil_dir_idx = findall(x -> x == true, isairfoil)

    old_path = path[dir_idx][airfoil_dir_idx]
    new_path = joinpath.(results_path, old_path)
    for (op, np) in zip(old_path, new_path)
        mv(op, np)
    end
end

function clear_directories(path, file_idx::Vector, dir_idx::Vector, field_file_idx::Vector, field_dir_idx::Vector, cellfields::Dict)
    for (i, field_name) in enumerate(keys(cellfields))
        if field_name == "uh"
            file_names = path[file_idx][field_file_idx[i]]
            dir_names = path[dir_idx][field_dir_idx[i]]
            map(x -> rm(x, recursive=true), file_names)
            map(x -> rm(x, recursive=true), dir_names)
        elseif field_name == "n_Γ"
            file_names = path[file_idx][field_file_idx[i]][2:end]
            dir_names = path[dir_idx][field_dir_idx[i]][2:end]
            map(x -> rm(x, recursive=true), file_names)
            map(x -> rm(x, recursive=true), dir_names)
        else

            file_names_clean = map(x -> x[9:end-4], path[file_idx][field_file_idx[i]])
            dir_names_clean = map(x -> x[9:end], path[dir_idx][field_dir_idx[i]])


            directory_already_analyzed = intersect(file_names_clean, dir_names_clean)
            directory_2_delete = map(x -> joinpath(results_path, x), directory_already_analyzed)
            println("Cleaning $(field_name)")
            map(x -> rm(x, recursive=true), directory_2_delete)
        end
    end
end


#Nodes analysis
function get_nodes(path::Vector{String})
    v = occursin.("nodes", path)
    idx_n = findall(x -> x == true, v)
    folder_nodes = path[idx_n][1]
    nodes_files = readdir(folder_nodes)
    df_nodes = DataFrame()
    for fnodes in nodes_files
        file_path = joinpath(folder_nodes, fnodes)
        df_tmp = DataFrame(CSV.File(file_path))
        df_nodes = vcat(df_nodes, df_tmp)
    end


    nodes_unique = unique(df_nodes)

    unique_idx = findall(Bool.(1 .- nonunique(df_nodes)))

    file_write = joinpath(results_path, "nodes_file.csv")

    CSV.write(file_write, df_nodes[unique_idx, :])

    return nodes_unique, unique_idx

end

#Export normals
function get_normals(path::Vector{String},unique_idx)
    v = occursin.("n_Γ", path)
    idx_n = findall(x -> x == true, v)
    folder_normals = path[idx_n][1]
    normals_files = readdir(folder_normals)
    df_normal = DataFrame()
    for fnormal in normals_files
        file_path = joinpath(folder_normals, fnormal)
        df_tmp = DataFrame(CSV.File(file_path))
        df_normal = vcat(df_normal, df_tmp)
    end

    normals_unique =  df_normal[unique_idx, :]
    file_write = joinpath(results_path, "normal_file.csv")

    CSV.write(file_write, normals_unique)

return normals_unique
end



function read_field_directories(path::Vector{String}, dir_idx::Vector, field_dir_idx::Vector, cellfields::Dict, unique_idx::Vector)
    for (field_number, field_name) in enumerate(keys(cellfields))
        field_directories = path[dir_idx][field_dir_idx[field_number]]

        for fd in field_directories[1:end]
            _, th = custom_cmp_dir(fd)
            file_write = joinpath(results_path, "$(field_name)_$(th).csv")
            field_files = readdir(fd)
            df_field = DataFrame()
            for f_file in field_files
                file_path = joinpath(fd, f_file)

                df_tmp = DataFrame(CSV.File(file_path))
                df_field = vcat(df_field, df_tmp)

            end
            println(file_write)

            CSV.write(file_write, df_field[unique_idx, :])

        end
    end
end


function average_field(path::Vector{String}, field_name::String, cellfields::Dict, file_idx::Vector, field_file_idx::Vector, unique_idx::Vector; offset=1, offend = 0)
    value_type, field_num = cellfields[field_name]
if offend>0
    vector_files = path[file_idx][field_file_idx[field_num]][offset:offend]
else
    vector_files = path[file_idx][field_file_idx[field_num]][offset:end]

end
    zv = zeros(length(unique_idx))

    if value_type == :scalar_value
        df_field = DataFrame(p=zv)

    elseif value_type == :vector_value
        df_field = DataFrame(x=zv, y=zv, z=zv)
    else
        error("type not recognized")
    end



    n_time_steps = length(vector_files)

    for file_read in vector_files
        df_tmp = DataFrame(CSV.File(file_read))
        println(file_read)
        df_field = df_field .+ df_tmp ./ n_time_steps
    end

    #Option for 3D average
    nodes_file_path = results_path * "nodes_file.csv"
    dfnodes = DataFrame(CSV.File(nodes_file_path))
    if !isempty(findall(z-> z > 0, dfnodes.z))
        println("3D Averaging")
        Zpoints = find_z_aligned(dfnodes)
        NZ = length(Zpoints[1])

        zv = zeros(length(Zpoints))
        df_field_avg =  DataFrame(p=zv)


        for (xyidx,zidx) in enumerate(Zpoints)
            z_avg = sum(df_field[zidx,:p]) / NZ
            df_field_avg[xyidx,:p] = z_avg
        end

    else
        df_field_avg = df_field

    end

    return df_field_avg
end


function extract_airfoil_features(nodes_unique::DataFrame, n_Γ0::DataFrame, Ph::DataFrame, Friction::DataFrame; u0::Float64, A::Float64, rho::Float64, α::Float64)
    q = 0.5 .* u0^2 * rho
    p1 = [0.0, 0.0] #leading edge
    p2 = [cosd(α), -1 * sind(α)] #trailing edge


 

    idx_nodes_airfoil = findall(abs.(nodes_unique.y) .< 0.5)
    nodes_airfoil0 = nodes_unique[idx_nodes_airfoil, :]
    n_Γ_airfoil0 = n_Γ0[idx_nodes_airfoil, :]

    if !isempty(findall(z-> z > 0, nodes_airfoil0.z))

    
    Zpoints = find_z_aligned(nodes_airfoil0)
    
    XYunique = map(x->x[1], Zpoints)

    
    nodes_airfoil = nodes_airfoil0[XYunique,:]
    n_Γ_airfoil = n_Γ_airfoil0[XYunique,:]
else
    nodes_airfoil = nodes_airfoil0
    n_Γ_airfoil = n_Γ_airfoil0
    end
    
    
    XY_Airfoil = Vector[]
    for (i, j) in zip(nodes_airfoil.x, nodes_airfoil.y)
        push!(XY_Airfoil, [i, j])
    end



    rev_idx = findall(is_above.(XY_Airfoil;p1,p2) .* sign.(n_Γ_airfoil.y) .<0) #reverse n_Γ sign for this

    n_Γ_airfoil.y[rev_idx] = -1 .* n_Γ_airfoil.y[rev_idx] 

    idx_top = findall(n_Γ_airfoil.y .> 0)
    idx_bottom = findall(n_Γ_airfoil.y .< 0)

    top_nodes = nodes_airfoil[idx_top, :]
    perm_top = sortperm(top_nodes.x)
    top_nodes_perm = top_nodes[perm_top, :]


    bottom_nodes = nodes_airfoil[idx_bottom, :]
    perm_bottom = sortperm(bottom_nodes.x)
    bottom_nodes_perm = bottom_nodes[perm_bottom, :]


    cp_top = Ph[idx_top, :][perm_top, :] ./ q
    cp_bottom = Ph[idx_bottom, :][perm_bottom, :] ./ q

    friction_top = .-Friction[idx_top, :][perm_top, :] ./ q .* μ
    friction_bottom = Friction[idx_bottom, :][perm_bottom, :] ./ q .* μ



    # cp_top = Ph[idx_nodes_airfoil, :][idx_top, :][perm_top, :] ./ q
    # cp_bottom = Ph[idx_nodes_airfoil, :][idx_bottom, :][perm_bottom, :] ./ q

    # friction_top = .-Friction[idx_nodes_airfoil, :][idx_top, :][perm_top, :] ./ q .* μ
    # friction_bottom = Friction[idx_nodes_airfoil, :][idx_bottom, :][perm_bottom, :] ./ q .* μ

    n_Γ_airfoil_top = n_Γ_airfoil[idx_top, :][perm_top, :]
    n_Γ_airfoil_bottom = n_Γ_airfoil[idx_bottom, :][perm_bottom, :]


    Float64.(top_nodes_perm.x), Float64.(bottom_nodes_perm.x),
    Float64.(top_nodes_perm.y), Float64.(bottom_nodes_perm.y),
    Float64.(cp_top.p), Float64.(cp_bottom.p), Float64.(friction_top.p), Float64.(friction_bottom.p), n_Γ_airfoil_top, n_Γ_airfoil_bottom

end





"""
It provides a Vector of length = airfoil points in 2D. Each element is a vector, where the elements have the same x and y.
"""
function find_z_aligned(df::DataFrame)
    zpoints = findall(x-> x == 0, df.z)

    
    ZPoints = Vector[]
    
    for zidx in zpoints
        xidx = findall(x->x== df[zidx,:].x, df.x)
        yidx = findall(y->y== df[zidx,:].y, df.y)
        line_points = intersect(xidx,yidx) 
        push!(ZPoints,line_points)
    end
    
    ZPoints
    
end

function get_tangent_x(V::Vector)
    x, y, z = V
    if y > 0
        t = [y, -x, z]
    else
        t = [-y, x, z]
    end
    return t
end


# function is_above(p; p1, p2)
#     m = (p2[2] - p1[2]) / (p2[1] - p1[1])
#     q = -m * p1[1] + p1[2]
#     res = p[2] - m * p[1] + q
#     if res > 0
#         return 1
#     else
#         return -1
#     end
# end

function is_above(p; p1, p2)
der = p2[2]/p2[1] .* 1.5#derivative in (0,0) #du89 AoA 1

c = -der/8 + der #derivative in (0,0) #du89 AoA 1

# c = -der/8 + der #derivative in (0,0)

a = (2*p2[2]-p2[1]*der - c*p2[1])/(-p2[1]^3)

b =(der-c-3*a*p2[1]^2)/(2 * p2[1])

treshold_fun(x) = a*x^3 +b*x^2 +c*x

    res = p[2] - treshold_fun(p[1])
    if res > 0
        return 1
    else
        return -1
    end
end

function read_cfd_xlsx(filename::String)
coeff_top = DataFrame(XLSX.readtable(filename,"top"))
coeff_bottom = DataFrame(XLSX.readtable(filename,"bottom"))
coeff_top[!,:] = convert.(Float64,coeff_top[!,:])
coeff_bottom[!,:] = convert.(Float64,coeff_bottom[!,:])

return coeff_top, coeff_bottom
end


### Not used

function compute_CDp()
    idx_inv_top = findfirst((top_nodesy[2:end] .- top_nodesy[1:end-1]) .< 0)
idx_inv_bottom = findfirst((bottom_nodesy[2:end] .- bottom_nodesy[1:end-1]) .> 0)

y_forw = [reverse(bottom_nodesy[1:idx_inv_bottom]);top_nodesy[1:idx_inv_top]]
y_back = [bottom_nodesy[idx_inv_bottom+1:end];reverse(top_nodesy[idx_inv_top+1:end])]
x_forw = [reverse(bottom_nodesx[1:idx_inv_bottom]);top_nodesx[1:idx_inv_top]]
x_back =[bottom_nodesx[idx_inv_bottom+1:end];reverse(top_nodesx[idx_inv_top+1:end])]

Cp_forw =  [reverse(cp_bottom[1:idx_inv_bottom]);cp_top[1:idx_inv_top]]
Cp_back = [cp_bottom[idx_inv_bottom+1:end];reverse(cp_top[idx_inv_top+1:end])]
return trapz(y_forw,Cp_forw)-trapz(y_back,Cp_back)
end