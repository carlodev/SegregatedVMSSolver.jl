module WallDistance

using Gridap
using LinearAlgebra: norm
using GridapGmsh
using GridapPETSc
using DataFrames
using CSV
using SegregatedVMSSolver.ExportUtility

import SegregatedVMSSolver.ExportUtility: conv_VectorValue


export get_initial_conditions

"""
    boundary_layer_height(Re::Real, L::Float64)

From Reynolds and characteristic dimension extimates the boundary layer height
"""
function boundary_layer_height(Re::Real, L::Float64)
  if Re < 500e3
    δ99 = 4.91 * L / (Re^0.5)

  else
    δ99 = 0.38 * L / (Re^0.2)

  end
  return δ99

end

"""
    boundary_layer_quadratic_profile(x, δ99::Float64, u_in::Real)

Parabolic approximation for boundary layer growth
"""
function boundary_layer_quadratic_profile(x, δ99::Float64, u_in::Real)
  if x > δ99
    return u_in
  else
    return (-(x / δ99) .^ 2 + 2 * (x / δ99)) * u_in
  end

end



"""
    intialize_picard_iteration(solver, Vg, V0, dΩ)

Solve the p poisson for p=1
"""
function intialize_picard_iteration(solver, Vg, V0, dΩ)
  flux(∇u) = ∇u
  f(x) = 1
  a(u, v) = ∫(∇(v) ⊙ (flux ∘ ∇(u))) * dΩ
  l(v) = ∫(v * f) * dΩ

  op = AffineFEOperator(a, l, Vg, V0)

  uh = interpolate_everywhere(1.0, Vg)

  uh, cache = GridapPETSc.solve!(uh, solver, op)

  return uh
end


"""
    run_picard_iteration(solver, Vg, V0, dΩ, uh_adv, p, γ)

Solve the linearized p-Poisson. `uh_adv` is the solution at the previous iteration. `γ` is the relaxation parameter.
"""
function run_picard_iteration(solver, Vg, V0, dΩ, uh_adv, p::Int64, γ::Float64)

  flux(∇u, ∇uha) = norm(∇uha)^(p - 2) * ∇u
  f(x) = 1
  a(u, v) = ∫(∇(v) ⊙ (flux ∘ (∇(u), ∇(uh_adv)))) * dΩ
  l(v) = ∫(v * f) * dΩ

  uh = interpolate_everywhere(1.0, Vg)

  op = AffineFEOperator(a, l, Vg, V0)
  uh, cache = GridapPETSc.solve!(uh, solver, op)


  return uh * γ + uh_adv * (1 - γ)
end



function initial_velocity_wall_distance((U0,Ug), (V0,Vg), Ω, dΩ, p::Int64, γ::Float64, u_in::Float64, δ99::Float64)

  options = "-ksp_type gmres -pc_type gamg -snes_rtol 1.0e-10 -ksp_converged_reason -ksp_error_if_not_converged true" #-ksp_monitor



  GridapPETSc.with(args=split(options)) do


    solver = PETScLinearSolver()

    uh2 = intialize_picard_iteration(solver, Vg, V0, dΩ)

    uh3 = run_picard_iteration(solver, Vg, V0, dΩ, uh2, 3, 1.0)
    uh = uh3

    for pi in 4:p
      for inner_iter in 1:1:3
        uh = run_picard_iteration(solver, Vg, V0, dΩ, uh, pi, γ)
        println("p = $pi, inner_iter = $(inner_iter)")
      end
    end

    vnorm(up, ∇up) = -norm(∇up)^(p - 1) + (p / (p - 1) * up + norm(∇up)^(p))^((p - 1) / p)
    vh = vnorm ∘ (uh, ∇(uh))

    DistFeFun = interpolate_everywhere(vh, Vg)

    function u_star_fun(d)
      xvel = boundary_layer_quadratic_profile(d, δ99, u_in)
      return VectorValue(xvel, 0.0)
    end


    u_start_oper = u_star_fun ∘ DistFeFun
    u_start_FEFun = interpolate_everywhere(u_start_oper, Ug)

    writevtk(Ω, "Inital_Condition", cellfields=["uh" => u_start_FEFun, "vh_norm" => vh, "uh_p6" => uh])

    u_start = get_free_dof_values(u_start_FEFun)

    return u_start, u_start_FEFun, Ω

  end #end GridapPETSc.do
end


function get_wall_distance(mesh_file::String,D::Int64, u_in::Float64, Re::Real,chord::Float64, walltag, p::Int64, γ::Float64)


  @assert D == 2 "Find wall distance not supported for 3D cases"
  @assert p>3 "p needs to be al least 4"
  @assert 0.0<γ<1.0 "Relaxation factor γ not valid"

  mesh_path = joinpath(mesh_file)
  model = GmshDiscreteModel(mesh_path)


  order = 1
  reffe = ReferenceFE(lagrangian, Float64, order)

  U0 = TestFESpace(model, ReferenceFE(lagrangian, VectorValue{D,Float64}, order); conformity=:H1, dirichlet_tags=walltag)
  Ug = TrialFESpace(U0, VectorValue(0.0, 0.0))

  V0 = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags=walltag)
  Vg = TrialFESpace(V0, zeros(length(walltag)))

  degree = 2
  Ω = Triangulation(model)
  dΩ = Measure(Ω, degree)

  δ99 = boundary_layer_height(Re, chord)

  
  u_start, u_start_oper, u_start_FEFun = initial_velocity_wall_distance((U0,Ug), (V0,Vg), Ω, dΩ, p, γ, u_in, δ99)

  return u_start, u_start_oper, u_start_FEFun

end



function export_intial_csv(Ω, uh_FE)
  f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)

  ref_grids = map(f, Gridap.Geometry.get_reffes(Ω))
  visgrid = Gridap.Visualization.VisualizationGrid(Ω, ref_grids)
  visgrid_ = conv_VectorValue.(visgrid.sub_grid.node_coordinates)
  unique_idx = unique(i -> visgrid_[i], eachindex(visgrid_))

  nodes_model = visgrid_[unique_idx] #Coordinate of unique node

  pdata = Gridap.Visualization._prepare_pdata(Ω, Dict("uh" => uh_FE), visgrid.cell_to_refpoints)
  uh_start = pdata["uh"][unique_idx]

  df_nodes = conv_to_df(nodes_model)
  df_uh = conv_to_df(uh_start)
  rename!(df_uh, ["uh_0", "uh_1", "uh_2"])

  df_export = hcat(df_nodes, df_uh)
  filename = joinpath("BoundaryLayerInit.csv")
  CSV.write(filename, df_export)
  return df_export
end

"""
  get_initial_conditions(mesh_file::String, u_in::Float64, Re::Real, walltag; chord::Float64=1.0, D::Int64=2, p::Int64=6)

It runs Picard iterations of the p-Poisson problem up to p=6 with relaxation factor of 0.5.
It is require p>3
It creates the `.vtu` file with the boundary layer initialization and also the `.csv` file which can be used to start a simulation over an airfoil.
It return a dataframe where for each nodes are provided the coordinates and the velocity vector components.

"""
function get_initial_conditions(mesh_file::String, u_in::Float64, Re::Real, walltag; chord::Float64=1.0, D::Int64=2, p::Int64=6, γ::Float64=0.5)
  u_start, uh_FE, Ω = get_wall_distance(mesh_file, D, u_in, Re, chord, walltag,p,γ)
  df_start = export_intial_csv(Ω, uh_FE)
  return df_start
end


end #end module