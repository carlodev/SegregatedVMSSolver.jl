using CSV, DataFrames, Plots, XLSX
using UnicodePlots
using Trapz
include("ReadResults_Utils.jl")

Re = 500_000
u0 = 45.0/2
c = 0.2
rho = 1.0
μ = u0*c*rho/Re
α = 1.0

cd("..")
cd("utils")
cd("SegregatedSolver")

cd("run_DU89_mesh_study/A1_500e3_C")

cd("run_DU89_mesh_study/A1_500e3_Fa")
cd("run_DU89_mesh_study/A1_500e3_Fz")
cd("run_DU89_mesh_study/A1_500e3_dt")
cd("run_DU89_mesh_study/A1_500e3_dt05")
cd("run_DU89_mesh_study/A1_500e3_dt2")
cd("run_DU89_mesh_study/A1_500e3_SSF")

cd("run_DU89_3D_A1_250e3")
cd("run_DU89_3D_A5_250e3")

cd("run_DU89_3D_A5_500e3")
cd("run_DU89_WT")
cd("run_DU89_WT_5")
cd("run_DU89_WT_250e3")
cd("A5")

pwd()

cd("run_params_study")

global results_path = "Results/"
# dir_results = ["T3"]

# for dr in dir_results

cellfields = Dict("ph"=>(:scalar_value,1), "friction"=>(:scalar_value,2))


file_idx, dir_idx, path = get_file_dir_idx(results_path)
field_file_idx, field_dir_idx = get_file_dir_idx_fields(path, file_idx, dir_idx, cellfields)

nodes_unique,unique_idx = get_nodes(path)
n_Γ = get_normals(path,unique_idx)



read_field_directories(path, dir_idx, field_dir_idx, cellfields, unique_idx)


file_idx, dir_idx, path = get_file_dir_idx(results_path)
field_file_idx, field_dir_idx = get_file_dir_idx_fields(path, file_idx, dir_idx, cellfields)

clear_directories(path,file_idx,dir_idx,field_file_idx, field_dir_idx,cellfields)


cd("..")

#ReadFiles
Ph = average_field(path, "ph", cellfields, file_idx, field_file_idx, unique_idx;offset = 5000, offend=0)
Friction = average_field(path, "friction", cellfields, file_idx, field_file_idx, unique_idx;offset = 5000, offend = 0)

include("ReadResults_Utils.jl")

top_nodesx,bottom_nodesx,top_nodesy,bottom_nodesy,cp_top,cp_bottom,friction_top,friction_bottom,n_Γ_airfoil_top,n_Γ_airfoil_bottom= extract_airfoil_features(nodes_unique, n_Γ, Ph, Friction; u0=u0, A=c, rho=rho, α = α)

t_Γ_airfoil_top = map(x->get_tangent_x([n_Γ_airfoil_top[x,:]...]), 1:1:size(n_Γ_airfoil_top)[1] )
t_Γ_airfoil_bottom = map(x->get_tangent_x([n_Γ_airfoil_bottom[x,:]...]), 1:1:size(n_Γ_airfoil_bottom)[1] )

scatter(top_nodesx,top_nodesy)
scatter!(bottom_nodesx,bottom_nodesy)


fname = "DU89_AoA$(Int64(α))_$(Int64(Re))"

starccm_file = "STARCCM+_"*fname*".xlsx"
# starccm_file = DataFrame(XLSX.read("Calderer.xlsx","top"))

# topr = DataFrame(XLSX.readtable("Calderer.xlsx", "top"))
# bottomr = DataFrame(XLSX.readtable("Calderer.xlsx", "bottom"))

coeff_top, coeff_bottom = read_cfd_xlsx(starccm_file)

plotly()
plt_Cf = plot(title=fname, xlabel ="x/c",ylabel="Cf",ylims=([-0.02,0.02]))
plot!(top_nodesx ./c, friction_top ,linecolor =:red, label = "VMS")
plot!(bottom_nodesx ./c, friction_bottom ,linecolor =:red,label = false)

# plot!(Float64.(topr.x), Float64.(topr.Cf),linecolor =:black,label = false)
# plot!(Float64.(bottomr.x), Float64.(bottomr.Cf),linecolor =:black,label = false)



plot!(coeff_top.xf,coeff_top.Cf, linecolor =:blue, label = "STARCCM+")
plot!(coeff_bottom.xf,coeff_bottom.Cf, linecolor =:blue, label = false)


Plots.savefig(plt_Cf, "Cf_"*fname*"_3D.pdf")


plotly()
plt_Cp = plot(title=fname,xlabel ="x/c",ylabel="Cp")
plot!(top_nodesx ./c, cp_top .+0.42,linecolor=:red, label = "VMS")
plot!(bottom_nodesx ./c, cp_bottom .+0.42,linecolor =:red, label = false)

yflip!()
# plot!(topr.xp, -topr.Cp,linecolor =:black,label = false)

plot!(coeff_top.xp ,coeff_top.Cp .+ 0.05 , linecolor =:blue, label = "STARCCM+")
plot!(coeff_bottom.xp,coeff_bottom.Cp .+ 0.05, linecolor =:blue, label = false)
yflip!()

Plots.savefig(plt_Cp, "Cp_"*fname*"_3D.pdf")


using Interpolations, BSplineKit, NumericalIntegration

pb_L = BSplineKit.interpolate(bottom_nodesx, cp_bottom, BSplineOrder(2))
pt_L = BSplineKit.interpolate(top_nodesx, cp_top, BSplineOrder(2))







Clfun(x) = (pb_L(x)-pt_L(x))
xtr = collect(0.0:0.001:1.0)
trapz(xtr, Clfun.(xtr))


CL_p = (trapz(bottom_nodesx,cp_bottom) - trapz(top_nodesx,cp_top)) ./c

CD_p = trapz(top_nodesy,cp_top)  - trapz(bottom_nodesy,cp_bottom) + trapz([bottom_nodesy[1],top_nodesy[1]],[cp_bottom[1],cp_top[1]])

CD_f = trapz(top_nodesx,friction_top) + trapz(bottom_nodesx,friction_bottom) 


CD = (CD_p+CD_f)/c


CL_p/CD


using JLD2

jldsave(fname*".jld2"; top_nodesx,bottom_nodesx,top_nodesy,bottom_nodesy,cp_top,cp_bottom,
    friction_top,friction_bottom,n_Γ_airfoil_top,n_Γ_airfoil_bottom,
    t_Γ_airfoil_top,t_Γ_airfoil_bottom)
