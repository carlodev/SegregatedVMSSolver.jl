module SolveProblem
using Gridap
using GridapDistributed
using GridapPETSc
using SparseArrays
using PartitionedArrays
using Parameters

using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using SegregatedVMSSolver.InitialConditions
using SegregatedVMSSolver.MatrixCreation
using SegregatedVMSSolver.VectorsOperations
using SegregatedVMSSolver.ExportUtility
using SegregatedVMSSolver.Interfaces

export solve_case



function initialize_solve(simcase::SimulationCase,params::Dict{Symbol,Any})
  @unpack trials,tests = params
  U,P = trials

  t0,dt,save_sim_dir = get_field(simcase,[:t0,:dt,:save_sim_dir])

  uh0, ph0 = create_initial_conditions(simcase,params)
  @info "Initial Conditions Created"

  matrices = initialize_matrices(trials,tests, t0+dt, uh0, params,simcase)
  vectors = initialize_vectors(matrices,uh0,ph0)

  initialize_export_nodes(params, simcase)

  mkpath(save_sim_dir)


  uh_avg = FEFunction(U(t0), vectors[2])
  set_zeros!(uh_avg.fields)
  ph_avg = FEFunction(P(t0),  vectors[1])
  set_zeros!(ph_avg.fields)
  return matrices, vectors, (uh_avg,ph_avg)
end




"""
  solve_case(params::Dict{Symbol,Any})

It solves iteratively the velocity and pressure system.
"""
function solve_case(params::Dict{Symbol,Any},simcase::SimulationCase)
  @unpack trials,tests = params
  U,P = trials

t0,dt,tF =  get_field(simcase.simulationp.timep,[:t0,:dt,:tF])
petsc_options,matrix_freq_update,M,a_err_threshold,θ =  get_field(simcase.simulationp.solverp,[:petsc_options,:matrix_freq_update,:M,:a_err_threshold,:θ])

time_step = collect(t0+dt:dt:tF)

init_values = initialize_solve(simcase,params)

matrices, vectors, (uh_avg,ph_avg) =  init_values


GridapPETSc.with(args=split(petsc_options)) do


Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App, 
  Mat_ML, Mat_inv_ML, Mat_S, Vec_Au, Vec_Ap = matrices

vec_pm,vec_um,vec_am,vec_sum_pm,Δa_star,Δpm1,Δa,b1,b2,ũ_vector = vectors

ns1 = create_PETSc_setup(Mat_ML,vel_kspsetup)
ns2 = create_PETSc_setup(Mat_S,pres_kspsetup)
uh_tn_updt = FEFunction(U(t0), vec_um)


for (ntime,tn) in enumerate(time_step)

      m = 0

      if mod(ntime,matrix_freq_update)==0
        println("update_matrices")
    
        @time begin
         Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App, 
          Mat_ML, Mat_inv_ML, Mat_S, Vec_Au, Vec_Ap = compute_matrices(trials, tests, tn+dt, uh_tn_updt, params,simcase)
         
          matrices = ( Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App, 
          Mat_ML, Mat_inv_ML, Mat_S, Vec_Au, Vec_Ap)

        end

        println("update numerical set up")
        @time begin
          numerical_setup!(ns1,Mat_ML)
          numerical_setup!(ns2,Mat_S)
        end

      end



      time_solve = @elapsed begin 


        vec_am .= pazeros(Mat_ML)
        vec_sum_pm .= pazeros(Mat_Aup)

        norm_Δa0 = 10
        norm_Δp0 = 10
        err_norm_Δa0 = 1
        err_norm_Δp0 = 1

        
      while (m<= M) && (err_norm_Δa0<a_err_threshold)

        Δpm1 .=  pazeros(Mat_S)
        Δa_star .= pazeros(Mat_ML)

        vectors = (vec_pm,vec_um,vec_am,vec_sum_pm,Δa_star,Δpm1,Δa,b1,b2,ũ_vector)
        println("norm um = $(norm(vec_um))")
        println("norm pm = $(norm(vec_pm))")
        println("norm da = $(norm(Δa))")
        println("norm am = $(norm(vec_am))")

        solve_velocity!(ns1,matrices,vectors,dt,θ)

        solve_pressure!(ns2,matrices,vectors,dt,θ)
      
        Δpm1 = GridapDistributed.change_ghost(Δpm1, Mat_Aup)
        println("norm pm1 = $(norm(Δpm1))")

        Δa .= Δa_star - θ .* Mat_inv_ML .* (Mat_Aup * Δpm1)

        vec_um .+=  dt * Δa
        vec_pm .+= Δpm1

        println("inner iter = $m")
        if m == 0
          vec_sum_pm .= Δpm1
          vec_am .= Δa
          norm_Δa0 = norm(Δa)
          norm_Δp0 = norm(Δpm1)
        else
          vec_sum_pm .+= Δpm1
          vec_am .+= Δa
        end
        
        err_norm_Δa0 = norm_Δa0/norm(Δa)
        err_norm_Δp0 = norm_Δp0/norm(Δpm1)

        evaluate_convergence(err_norm_Δa0, "velocity")
        evaluate_convergence(err_norm_Δp0, "pressure")

        m = m + 1

    end  #end while
  end #end elapsed


  println("solution time")
  println(time_solve)
    GridapPETSc.GridapPETSc.gridap_petsc_gc()

update_ũ_vector!(ũ_vector,vec_um)

uh_tn_updt = FEFunction(U(tn+dt), vec_um)

if ntime>100
  uh_tn_updt = FEFunction(U(tn+dt), update_ũ(ũ_vector))
end

println("Solution computed at time $tn")
uh_tn = FEFunction(U(tn), vec_um)
ph_tn = FEFunction(P(tn), vec_pm)

update_time_average!(uh_tn,ph_tn, uh_avg, ph_avg, tn, time_step, simcase.simulationp.timep)

writesolution(params, simcase, ntime, tn, (uh_tn,ph_tn,uh_tn_updt,uh_avg,ph_avg))

export_fields(params,simcase, tn, uh_tn, ph_tn)


  end #end for
end #end GridapPETSc

end #end solve_case


function solve_velocity!(ns1, matrices, vectors, dt::Float64, θ::Float64)
  Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App, Mat_ML, Mat_inv_ML, Mat_S, Vec_Au, Vec_Ap = matrices
  vec_pm,vec_um,vec_am,vec_sum_pm,Δa_star,Δpm1,Δa,b1,b2,ũ_vector = vectors

  vec_um = GridapDistributed.change_ghost(vec_um, Mat_Auu)
  vec_pm = GridapDistributed.change_ghost(vec_pm, Mat_Aup)
  vec_am = GridapDistributed.change_ghost(vec_am, Mat_ML)

  println("solving velocity")
    
    b1 .= -Mat_Auu * vec_um - Mat_Aup * vec_pm - Mat_ML * vec_am +
    Mat_Auu * dt * vec_am + (1 - θ) * Mat_Aup * vec_sum_pm + Vec_Au

  @time solve!(Δa_star,ns1,b1)
end


function solve_pressure!(ns2, matrices, vectors, dt::Float64, θ::Float64)
  Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App, Mat_ML, Mat_inv_ML, Mat_S, Vec_Au, Vec_Ap = matrices
  vec_pm,vec_um,vec_am,vec_sum_pm,Δa_star,Δpm1,Δa,b1,b2,ũ_vector = vectors

  vec_um = GridapDistributed.change_ghost(vec_um, Mat_Apu)
  vec_pm = GridapDistributed.change_ghost(vec_pm, Mat_App)
  vec_am = GridapDistributed.change_ghost(vec_am, Mat_Tpu)
  Δa_star = GridapDistributed.change_ghost(Δa_star, Mat_Tpu)

  println("solving pressure")

  #-Vec_A because changing sign in the continuity equations
  b2 .= Mat_Tpu * Δa_star + Mat_Apu * (vec_um + dt * Δa_star) + Mat_App * vec_pm + Mat_Tpu * vec_am - Vec_Ap

  @time solve!(Δpm1,ns2,b2)

end


end #end module 