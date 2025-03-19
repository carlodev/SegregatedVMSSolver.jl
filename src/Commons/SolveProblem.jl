module SolveProblem
using Gridap
using GridapDistributed
using GridapPETSc
using SparseArrays
using PartitionedArrays
using Parameters
using Gridap.FESpaces

using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using SegregatedVMSSolver.CreateProblem
using SegregatedVMSSolver.MatrixCreation
using SegregatedVMSSolver.VectorsOperations
using SegregatedVMSSolver.ExportUtility
using SegregatedVMSSolver.Interfaces

export solve_case



function initialize_solve(simcase::SimulationCase,params::Dict{Symbol,Any})
  @unpack trials,tests = params
  U,P = trials
  
  @sunpack t0,dt,save_sim_dir = simcase


  Ut0 = U(t0)
  Pt0 = P(t0)

  Ut0_1 = U(t0+dt)
  Pt0_1 = P(t0+dt)

  merge!(params, Dict(:Utn => Ut0, :Ptn => Pt0, :Utn1 => Ut0_1, :Ptn1 => Pt0_1))


  uh0, ph0 = create_initial_conditions(simcase,params)

  @info "Initial Conditions Created"

  matrices = initialize_matrices(uh0, params,simcase)
  vectors = initialize_vectors(matrices,uh0,ph0)


  initialize_export_nodes(params, simcase)

  mkpath(save_sim_dir)


  uh_avg = FEFunction(Ut0, vectors[2])
  set_zeros!(uh_avg.fields)
  ph_avg = FEFunction(Pt0,  vectors[1])
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
  
  @unpack trials,tests = params
  
@sunpack t0,dt,tF =  simcase.simulationp.timep
@sunpack petsc_options,matrix_freq_update,M,a_err_threshold,θ,Number_Skip_Expansion = simcase.simulationp.solverp
time_step = collect(t0+dt:dt:tF)

init_values = initialize_solve(simcase,params)


@unpack Utn, Utn1, Ptn, Ptn1 = params

matrices, vectors, (uh_avg,ph_avg) =  init_values


GridapPETSc.with(args=split(petsc_options)) do


Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App, 
  Mat_ML, Mat_inv_ML, Mat_S, Vec_Auu, Vec_Aup, Vec_Apu, Vec_App, Vec_Au, Vec_Ap = matrices

vec_pm,vec_um,vec_am,vec_sum_pm,Δa_star,Δpm1,Δa,b1,b2,ũ_vector = vectors

ns1 = create_PETSc_setup(Mat_ML,vel_kspsetup)
ns2 = create_PETSc_setup(Mat_S,pres_kspsetup)

uh_tn_updt = FEFunction(Utn, vec_um)

for (ntime,tn) in enumerate(time_step)

      m = 0
     
      @info "inner iteration $m | outer iteration $ntime"

      if mod(ntime,matrix_freq_update)==0

        @info "updating matrix and vectors"
        @time update_all_matrices_vectors!(matrices, uh_tn_updt, params,simcase)
        @info "matrix and vectors updated"


        @info "update numerical set up"
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

        solve_velocity!(ns1,matrices,vectors,dt,θ)

        solve_pressure!(ns2,matrices,vectors,dt,θ)
      
     
        Δpm1 = GridapDistributed.change_ghost(Δpm1, Mat_Aup)

       Δa .= Δa_star - θ .* Mat_inv_ML .* (Mat_Aup * Δpm1)

        vec_um .+=  dt * Δa
        vec_pm .+= Δpm1

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


  println("solution time at time $tn")
  println(time_solve)
    @time GridapPETSc.GridapPETSc.gridap_petsc_gc()

update_ũ_vector!(ũ_vector,vec_um)



Utn = Utn1
Utn1 = U(tn+dt)

Ptn = Ptn1
Ptn1 = P(tn+dt)

params[:Ptn1] = Ptn1
params[:Utn1] = Utn1

uh_tn_updt = FEFunction(Utn1, vec_um)

if ntime>Number_Skip_Expansion
  uh_tn_updt = FEFunction(Utn1, update_ũ(ũ_vector))
end

uh_tn = FEFunction(Utn, vec_um)
ph_tn = FEFunction(Ptn, vec_pm)

uh_avg = update_time_average(uh_tn, uh_avg, Utn, tn, ntime, time_step, simcase.simulationp.timep)
ph_avg = update_time_average(ph_tn, ph_avg, Ptn, tn, ntime, time_step, simcase.simulationp.timep)



writesolution(params, simcase, ntime, tn, (uh_tn,ph_tn), (uh_avg,ph_avg))

export_fields(params,simcase, tn, uh_tn, ph_tn)

  end #end for
end #end GridapPETSc

end #end solve_case


function solve_velocity!(ns1, matrices, vectors, dt::Float64, θ::Float64)
  Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App, Mat_ML, Mat_inv_ML, Mat_S,Vec_Auu, Vec_Aup, Vec_Apu, Vec_App, Vec_Au, Vec_Ap = matrices
  vec_pm,vec_um,vec_am,vec_sum_pm,Δa_star,Δpm1,Δa,b1,b2,ũ_vector = vectors

  vec_um = GridapDistributed.change_ghost(vec_um, Mat_Auu)
  vec_pm = GridapDistributed.change_ghost(vec_pm, Mat_Aup)
  vec_am = GridapDistributed.change_ghost(vec_am, Mat_ML)

  b1 .= -Mat_Auu * vec_um - Mat_Aup * vec_pm - Mat_ML * vec_am +
    Mat_Auu * dt * vec_am + (1 - θ) * Mat_Aup * vec_sum_pm + Vec_Au
    println("solving velocity")
  @time solve!(Δa_star,ns1,b1)

end


function solve_pressure!(ns2, matrices, vectors, dt::Float64, θ::Float64)
  Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App, Mat_ML, Mat_inv_ML, Mat_S, Vec_Auu, Vec_Aup, Vec_Apu, Vec_App, Vec_Au, Vec_Ap = matrices
  vec_pm,vec_um,vec_am,vec_sum_pm,Δa_star,Δpm1,Δa,b1,b2,ũ_vector = vectors

  vec_um = GridapDistributed.change_ghost(vec_um, Mat_Apu)
  vec_pm = GridapDistributed.change_ghost(vec_pm, Mat_App)
  vec_am = GridapDistributed.change_ghost(vec_am, Mat_Tpu)
  Δa_star = GridapDistributed.change_ghost(Δa_star, Mat_Tpu)


  #-Vec_A because changing sign in the continuity equations
  b2 .= Mat_Tpu * Δa_star + Mat_Apu * (vec_um + dt * Δa_star) + Mat_App * vec_pm + Mat_Tpu * vec_am - Vec_Ap
  println("solving pressure")
  @time solve!(Δpm1,ns2,b2)

end


end #end module 