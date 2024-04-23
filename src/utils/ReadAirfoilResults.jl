module ReadAirfoilResults

using CSV
using DataFrames
using Trapz
using Statistics
using FFTW
using ScatteredInterpolation



export get_nodes
export get_normals
export average_field
export extract_airfoil_features
export compute_CL_CD

export compute_average
export compute_average_2D
export compute_plane_tke
export compute_scatter_interp



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







###############################################
## Read Fluctuations
###############################################
"""
    ZProbe
Give a 2D point coordinate xp, yp, it gived the idx of the points aligned in the z direction
"""
struct ZProbe
    xp::Float64
    yp::Float64
    idx::Vector{Int64}
end

function ZProbe(df::DataFrame, xpyp::Vector{Float64})
    xp,yp=xpyp
    
    df0 =filter(x->x.z ==0, df)
    
    dx = df0.x .- xp 
    dy = df0.y .- yp
    
    _,idx0 = findmin(dx.^2 .+dy.^2)
    xc0 = df0[idx0,:].x
    yc0 = df0[idx0,:].y
    

    dfM=Matrix(df)

    idx0x = findall(dfM[:,1].== xc0)
    idx0y = findall(dfM[:,2].== yc0)
    @assert idx0x==idx0y
    
   
    z = dfM[idx0x,3]
    zsort = sortperm(z)
    ZProbe(xc0,yc0,idx0x[zsort])
end

"""
    XYPlane
Give a zp coordinate, it gives the indexes of the points in that XY plane
"""
struct XYPlane
    zp::Float64
    idx::Vector{Int64}
end

function XYPlane(df::DataFrame, zp::Float64)
    zidx = findall(x-> x == zp, df.z)
    tsort = zidx[sortperm(df[zidx,:])]
    # tsort_filter = tsort[df[tsort,:].x .>0.5]
    XYPlane(zp,tsort)
end


"""
    get_idx_sort_reduct(res_path::String,tagname::String,offset::Int64,offend::Int64)

In res_path directory, for the boundary tagname, it provides the sorted file indexes
"""
function get_idx_sort_reduct(res_path::String,tagname::String,offset::Int64,offend::Int64)
    filenames = readdir(res_path)
    idx_tag = findall(occursin.("$(tagname)_uh",filenames))

    idx_sort = sortperm(filenames[idx_tag], by=custom_cmp_file)

    if offend<0
        offend = length(idx_sort)
    else
        @assert offend>offset
        @assert offend<=length(idx_sort)
    end
    idx_sort_reduct = idx_tag[idx_sort[offset:offend]]

    return filenames, idx_sort_reduct
end


"""
    read_fluctuations(res_path::String, xpyp::Vector{Float64};tagname="topairfoil", offset=1,offend=-1)

It is reading the results in `res_path`, for all the points aligned in Z direction of coordinates `xpyp`. 
It provides the ZProbe which has the information of the actual nodes used.
Vel_Mat is a dense matrix: Time x Points x Velocity Components (3)
offset and offend can be used to skip intial and final files (avoiding OutOfMemory() error)
"""
function read_fluctuations(res_path::String, xpyp::Vector{Float64};tagname="topairfoil", offset=1,offend=-1)
    
    filenames, idx_sort_reduct = get_idx_sort_reduct(res_path,tagname,offset,offend)
    tnodes = get_nodes(res_path; tagname=tagname)

    zprobe = ZProbe(tnodes,xpyp)

    println("Probe in point $(zprobe.xp), $(zprobe.yp), idxs=$(zprobe.idx) ")
   
    Vel_Mat = zeros(length(idx_sort_reduct),length(zprobe.idx),3) #time x point x velocity direction


    for (i,fi) in enumerate(idx_sort_reduct)
        fpath = joinpath(res_path,filenames[fi])
        println(fpath)
        df_tmp = CSV.read(fpath, DataFrame)[zprobe.idx,:]
        Vel_Mat[i,:,:]=copy(Matrix(df_tmp))
        
    end
    

return zprobe,Vel_Mat
end

"""
    compute_PSD(Vel::Array, tn::Vector{Float64})

It computes the PSD of the TKE. Vel is the result obtained through `read_fluctuations`. The spectras are averaged in the Z direction 
"""
function compute_PSD(Vel::Array, tn::Vector{Float64})
    VelMean = Statistics.mean(Vel,dims=1)

    VelFluct = (Vel .- VelMean).^2

    Ktime = 0.5 .* sum(VelFluct,dims=3)
    Nz = size(Ktime)[2]

    freqs_z = Vector[]
    PSD_z = Vector[]

    dt = tn[2] - tn[1]
    nt=length(tn)


    for i = 1:1:Nz
        fhat=fft(Ktime[:,i,1])
        PSD = fhat.*conj(fhat)/(nt)
        PSD = real(fftshift(PSD))
        freqs = fftshift(fftfreq(nt,1/dt))
        push!(PSD_z,PSD)
        push!(freqs_z,freqs)

    end

    return mean(PSD_z), freqs_z[1]
end



"""
    compute_plane_tke(res_path::String;tagname="topairfoil", offset=1,offend=-1, zp=[0.1])

It computes the TKE for each point in the plane with Z=Zp.
"""
function compute_plane_tke(res_path::String;tagname="topairfoil", offset=1,offend=-1, zp=[0.1])


    filenames, idx_sort_reduct = get_idx_sort_reduct(res_path,tagname,offset,offend)
    tnodes = get_nodes(res_path; tagname=tagname)
    
    Zidx = Vector[]
    
    nz = length(zp)
    
    for z in zp
        zidx = XYPlane(tnodes, z).idx
        push!(Zidx,zidx)
    end
    
        
    Vel_Mat = zeros(nz,length(idx_sort_reduct), length(Zidx[1]),3) #number of planes x times  x pointx x velocity direction
    
    for (i,fi) in enumerate(idx_sort_reduct)
        fpath = joinpath(res_path,filenames[fi])
        println(fpath)
        for (j,zidx) in enumerate(Zidx)
            df_tmp = CSV.read(fpath, DataFrame)[zidx,:]
            Vel_Mat[j,i,:,:]=copy(Matrix(df_tmp))
        end

 
        
    end
    
    Vel_uprime = Statistics.std(Vel_Mat,dims=2)
    Vel_uprime2 = Vel_uprime .^2 
    Vel_uprime2_zavg = Statistics.mean(Vel_uprime2,dims=1)
    TKE = vec(0.5 .*sum(Vel_uprime2_zavg,dims=4))

    return TKE, Vel_Mat
    
end

##### Time Average

function compute_average(res_path::String;tagname="topairfoil", offset=1,offend=-1)
    
    filenames, idx_sort_reduct = get_idx_sort_reduct(res_path,tagname,offset,offend)
    Ntime = length(idx_sort_reduct)

    fpath = joinpath(res_path,filenames[idx_sort_reduct[1]])
    df_avg = CSV.read(fpath, DataFrame)./Ntime

    for (i,fi) in enumerate(idx_sort_reduct[2:end])
        fpath = joinpath(res_path,filenames[fi])
        println(fpath)
        df_tmp = CSV.read(fpath, DataFrame)
        df_avg = df_avg .+ df_tmp./Ntime
        
    end
    Vel_avg_3D=Matrix(df_avg)


    return Vel_avg_3D
end


function compute_average_2D(Vel_avg3D::Matrix; tagname="topairfoil")
    tnodes = get_nodes(res_path; tagname=tagname)

    unique_z = unique(tnodes.z)
    z_idx_sort = sortperm(unique_z)
    unique_z[z_idx_sort]
    
    zpoints = findall(x-> x == unique_z[z_idx_sort[2]], tnodes.z)
    Nzpoints = length(zpoints)
    Nz = length(z_idx_sort) - 2
    Vel_avg2D = zeros(Nzpoints,3)
    

    
    for z in  unique_z[z_idx_sort[2:end-1]]
        zidx = findall(x-> x == z, tnodes.z)
        tsort = sortperm(tnodes[zidx,:])

        M =     Vel_avg3D[zidx[tsort],:]
        Vel_avg2D = Vel_avg2D .+ M ./ Nz
    end

    return Vel_avg2D
end





####Interpolate


"""
    compute_scatter_interp(res_path, velocity::Vector{Float64} ,zp::Float64; tagname="topairfoil", ylims=[-0.018,0.15], xlims=[0.5,1.0])

In res_path is reading the nodes file associated with the tagname. It computes the interpolation of the values in velocity on nodes of plane in zp.
NearestNeighbor() algorithm is used.
"""
function compute_scatter_interp(res_path, velocity::Vector{Float64} ,zp::Float64; tagname="topairfoil", ylims=[-0.018,0.15], xlims=[0.5,1.0])

    tnodes = get_nodes(res_path; tagname=tagname)
    

    zidx = XYPlane(tnodes, zp).idx

    tz0nodes = tnodes[zidx,1:2]
    tz0nodesM=Matrix(tz0nodes)


    
    x = tz0nodes.x
    y = tz0nodes.y
    points = tz0nodesM'
    
    
    
    
    # Create dense grid of points
    x_dense = collect(range(xlims..., length=200))
    
    y_dense = vcat(collect((range(ylims..., length=250))))
    
    grid = collect(Iterators.product(x_dense, y_dense))
    xgrid = vec(getindex.(grid,1))
    ygrid = vec(getindex.(grid,2))
    
    
    itp = ScatteredInterpolation.interpolate(NearestNeighbor(), points, velocity)
    
    velocity_dense = zeros(length(grid))
    
    for (i,p) in enumerate(grid)
        pp = [p[1];p[2]]
        interpolated = ScatteredInterpolation.evaluate(itp, pp)
        velocity_dense[i] = interpolated[1]
    end
    
    return  xgrid,ygrid,velocity_dense
    
end




end #end module