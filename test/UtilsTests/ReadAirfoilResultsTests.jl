module ReadAirfoilResultsTests

using Test
using SegregatedVMSSolver
using SegregatedVMSSolver.ReadAirfoilResults
# using Plots

using DataFrames

res_path = joinpath(@__DIR__, "Results")
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

rm(res_path; force=true, recursive=true)

end #end module