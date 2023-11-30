module ReadAirfoilResultsTests

using Test
using SegregatedVMSSolver
using SegregatedVMSSolver.ReadAirfoilResults

using DataFrames

res_path = "Results/"


nodes = get_nodes(res_path; tagname="V75")
@test typeof(nodes) <: DataFrame

normals = get_normals(res_path)
@test typeof(normals) <: DataFrame

Ph = average_field(res_path, "ph", nodes)
Friction = average_field(res_path, "friction", nodes)

@test typeof(Ph) <: DataFrame
@test typeof(Friction) <: DataFrame

# plot(nodes.x,nodes.y,seriestype=:scatter)

Re = 500_000
u0 = 1.0
c = 1.0
rho = 1.0
μ = u0*c*rho/Re
α = 0.0


top_nodesx,bottom_nodesx,top_nodesy,bottom_nodesy,cp_top,cp_bottom,
friction_top,friction_bottom,n_Γ_airfoil_top,n_Γ_airfoil_bottom=extract_airfoil_features(nodes, normals, Ph, Friction; u0=u0, μ=μ, rho=rho, α=α, chord=c)

# plot(top_nodesx,top_nodesy, seriestype=:scatter)
# plot!(bottom_nodesx,bottom_nodesy, seriestype=:scatter)


CL,CD = compute_CL_CD(top_nodesx,bottom_nodesx,top_nodesy,bottom_nodesy,cp_top,cp_bottom,
friction_top,friction_bottom; chord = c)

@test typeof(CL)<:Real
@test typeof(CD)<:Real

rm(res_path; force=true, recursive=true)
end #end module