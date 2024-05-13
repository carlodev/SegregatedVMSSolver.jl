module ReadAirfoilResultsTests

using Test
using SegregatedVMSSolver
using SegregatedVMSSolver.ReadAirfoilResults
# using Plots

using DataFrames

for res_path_dir in ["Results2D","Results3D"]

    res_path = joinpath(@__DIR__, res_path_dir)
    Re = 500_000
    u0 = 1.0
    c = 1.0
    rho = 1.0
    μ = u0*c*rho/Re
    α = 1.0


    nodes, normals = get_geometry_info(res_path;α=α)

    
    @test typeof(nodes) <: GeometryNodes
    @test typeof(normals) <: GeometryNormals

    Ph = time_space_average_field(res_path, "ph", nodes)
    Velocity = time_space_average_field(res_path, "uh", nodes)
    Friction = time_space_average_field(res_path, "friction", nodes)
    
    @test typeof(Ph) <: DataFrame
    @test typeof(Friction) <: DataFrame
    @test typeof(Velocity) <: DataFrame

    Spanwise_Pressure_Avg = average_3D_field(res_path, "ph",nodes)
    Spanwise_Velocity_Avg = average_3D_field(res_path, "uh",nodes)
    Spanwise_Friction_Avg = average_3D_field(res_path, "friction",nodes)

        
    @test typeof(Spanwise_Pressure_Avg) <: Vector{DataFrame}
    @test typeof(Spanwise_Velocity_Avg) <: Vector{DataFrame}
    @test typeof(Spanwise_Friction_Avg) <: Vector{DataFrame}
    
    cp_top, cp_bottom = extract_Cp(nodes, Ph; u0=u0, rho=rho)
    friction_top, friction_bottom = extract_Cf(nodes, Friction, μ; u0=u0, rho=rho)

    @test typeof(cp_top) <: Vector{Float64}
    @test typeof(cp_bottom) <: Vector{Float64}
    @test typeof(friction_top) <: Vector{Float64}
    @test typeof(friction_bottom) <: Vector{Float64}

    CL,CD = compute_CL_CD(nodes, cp_top, cp_bottom,
    friction_top, friction_bottom; chord=1.0)


    @test typeof(CL)<:Real
    @test typeof(CD)<:Real

end #end for

###############################################
## Read Fluctuations Tests
###############################################
res_path_dir = "Results3D"
res_path = joinpath(@__DIR__, res_path_dir)

tagname="topairfoil"

topairfoil_nodes = get_nodes(res_path; tagname=tagname)

@test typeof(topairfoil_nodes) <: DataFrame

Vel_avg3D = compute_time_average(res_path; tagname=tagname)
Vel_avg2D = compute_time_span_average(res_path,Vel_avg3D; tagname=tagname)

@test typeof(Vel_avg3D) <:Matrix
@test typeof(Vel_avg2D) <:Matrix


ux_avg = Vel_avg2D[:,1]
uy_avg = Vel_avg2D[:,2]
uz_avg = Vel_avg2D[:,3]

@test typeof(ux_avg) <:Vector
@test typeof(uy_avg) <:Vector
@test typeof(uz_avg) <:Vector


xgrid,ygrid,velocity_dense = compute_scatter_interp(res_path, ux_avg, 0.1)

@test typeof(xgrid) <:Vector
@test typeof(ygrid) <:Vector
@test typeof(velocity_dense) <:Vector

TKE, Vel_Mat = compute_plane_tke(res_path; tagname=tagname)

@test typeof(TKE) <:Vector
@test typeof(Vel_Mat) <:Array

xc,yc= 0.5,0.1 #it has to be a point of the domain
zprobe, Vel = read_fluctuations(res_path, [xc,yc])

dt = 0.001
nt=size(Vel)[1]-1 #Number of time-steps
tn = collect(0.0:dt:dt*nt)
PSD,freqs = compute_PSD(Vel, tn)

@test typeof(PSD) <: Vector{Float64}
@test typeof(freqs) <: Vector{Float64}
@test length(PSD) == length(freqs)


end #end module