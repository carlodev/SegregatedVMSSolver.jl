module ReadAirfoilResults

using CSV
using DataFrames
using Trapz

export get_nodes
export get_normals
export average_field
export extract_airfoil_features
export compute_CL_CD



function custom_cmp_file(x::String)
    offset = findlast("_",x)[1]+1
    offend = findlast(".",x)[1]-1
    str, num_str = SubString(x, 1, offset-1), SubString(x, offset, offend)
    num = parse(Float64, num_str)
    return str, num
end

"""
    get_nodes(path::String)

It provides a `DataFrame` with the nodes of the `Airfoil` boundary
"""
function get_nodes(path::String; tagname="airfoil")
    f_path = readdir(path)
    idx_n = findfirst(x->occursin("$(tagname)_nodes", x), f_path)
    nodes_file_path = joinpath(path, f_path[idx_n])
    df_nodes = DataFrame(CSV.File(nodes_file_path))
    return df_nodes
end

"""
    get_nodes(path::String)

It provides a `DataFrame` with the normals vectors at the surface of the `Airfoil` boundary    
"""
function get_normals(path::String; tagname="airfoil")
    f_path = readdir(path)
    idx_n = findfirst(x->occursin("$(tagname)_n_Γ", x), f_path)
    normals_file_path = joinpath(path, f_path[idx_n])
    df_normals = DataFrame(CSV.File(normals_file_path))
    return df_normals
end


"""
    average_field(path::String, field_name::String, nodes::DataFrame; offset=1, offend_ = 0)

It computes the time-average and also the spanwise averge (z direction). `field_name` can be `ph` for pressure or `friction` for friction.
It is possible to skip a certain amount of initial time-step when averaging setting the `offset`.
"""
function average_field(path::String, field_name::String, nodes::DataFrame; offset=1, offend_ = 0, tagname="airfoil")
    f_path = readdir(path)
    file_name = tagname * "_" * field_name
    idx_n = findall(x->occursin(file_name, x), f_path)
    field_file_path = joinpath.(path, f_path[idx_n])
    
    offend = (offend_==0) ? length(field_file_path) : offend_

    @assert offset >0
    @assert offend > offset
    
    idx_sort = sortperm(field_file_path, by=custom_cmp_file)

    vector_files = field_file_path[idx_sort][offset:offend]
    n_time_steps = length(vector_files)
    
    df_field = DataFrame(CSV.File(vector_files[1])) ./ n_time_steps

    for fread in vector_files[2:end]
        println(fread)
        df_tmp = DataFrame(CSV.File(fread))
        df_field = df_field .+ df_tmp ./ n_time_steps
    end
   
    

    #3D averaging
    if !isempty(findall(z-> z > 0, nodes.z)) # In 2D case all the nodes have 0.0 as Z coordinates
        println("3D Averaging")
        Zpoints = find_z_aligned(nodes)
        NZ = length(Zpoints[1])

        # zv = zeros(length(Zpoints))
        df_field_avg =  copy(df_field)[1:length(Zpoints),:]


        for (xyidx,zidx) in enumerate(Zpoints)
            for pname in propertynames(df_field)
                z_avg = sum(df_field[zidx,pname]) / NZ
                df_field_avg[xyidx,pname] = z_avg
            end
        end

    else
        df_field_avg = df_field

    end

    return df_field_avg
end



# find_z_aligned(df::DataFrame)
# It provides a Vector of length = airfoil points in 2D. Each element is a vector, where the elements have the same x and y.
# 
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


"""
    extract_airfoil_features(nodes::DataFrame, n_Γ0::DataFrame, Ph::DataFrame, Friction::DataFrame; u0::Float64, μ::Float64, rho::Float64, α::Float64, chord::Float64)

It provides the airfoil features like pressure and friction coefficient splitted between top and bottom. The nodes and results are orderded for growing x coordinates.
"""
function extract_airfoil_features(nodes::DataFrame, n_Γ0::DataFrame, Ph::DataFrame, Friction::DataFrame; u0::Float64, μ::Float64, rho::Float64, α::Float64, chord::Float64, der_slope=-1.0)
    q = 0.5 .* u0^2 * rho
    p1 = [0.0, 0.0] #leading edge
    p2 = [chord *cosd(α), -chord * sind(α)] #trailing edge


 

    idx_nodes_airfoil = findall(abs.(nodes.y) .< 0.5)
    nodes_airfoil0 = nodes[idx_nodes_airfoil, :]
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



    rev_idx = findall(is_above.(XY_Airfoil;p1,p2, der_slope=der_slope) .* sign.(n_Γ_airfoil.y) .<0) #reverse n_Γ sign for this

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


    n_Γ_airfoil_top = n_Γ_airfoil[idx_top, :][perm_top, :]
    n_Γ_airfoil_bottom = n_Γ_airfoil[idx_bottom, :][perm_bottom, :]


    Float64.(top_nodes_perm.x), Float64.(bottom_nodes_perm.x),
    Float64.(top_nodes_perm.y), Float64.(bottom_nodes_perm.y),
    Float64.(cp_top.p), Float64.(cp_bottom.p), Float64.(friction_top.p), Float64.(friction_bottom.p), n_Γ_airfoil_top, n_Γ_airfoil_bottom

end

#  is_above(p; p1, p2)
#  It recognize if a specific point coordinate `p` is on the top or bottom side. It use a threshold line which is a cuibic function passing by the trailing `p2` and leading ´p1´ edges.
#  
function is_above(p; p1, p2, der_slope=-1)
    der = p2[2]/p2[1] .* -der_slope # -0.20 #derivative in trailing edge Manage the (-0.20) factor
    c = -der/8 + der #derivative in leading edge
    
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

"""
    compute_CL_CD(top_nodesx,bottom_nodesx,top_nodesy,bottom_nodesy,cp_top,cp_bottom,
    friction_top,friction_bottom; chord = 1.0)

It computes lift and drag coefficients
"""
function compute_CL_CD(top_nodesx,bottom_nodesx,top_nodesy,bottom_nodesy,cp_top,cp_bottom,
    friction_top,friction_bottom; chord = 1.0)

    CL = (trapz(bottom_nodesx,cp_bottom) - trapz(top_nodesx,cp_top)) ./chord

    #CD pressure
    CD_p = trapz(top_nodesy,cp_top)  - trapz(bottom_nodesy,cp_bottom) + trapz([bottom_nodesy[1],top_nodesy[1]],[cp_bottom[1],cp_top[1]])
    
    #CD friction
    CD_f = trapz(top_nodesx,friction_top) + trapz(bottom_nodesx,friction_bottom) 
    
    
    CD = (CD_p+CD_f)/chord
    
    return CL,CD
end


end #end module