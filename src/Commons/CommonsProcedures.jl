"""
  print_model(params::Dict{Symbol, Any})

It prints the model mesh
"""
function print_model(params::Dict{Symbol, Any})
  @unpack printmodel, model, case = params
  if printmodel
      writevtk(model,case)
  end
end


"""
  creation_fe_spaces(params::Dict{Symbol,Any}, u_diri_tags, u_diri_values, p_diri_tags, p_diri_values)

It creates the finite elements spaces accordingly to the previously generated dirichelet tags
"""
function creation_fe_spaces(params::Dict{Symbol,Any}, u_diri_tags, u_diri_values, p_diri_tags, p_diri_values)
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{params[:D],Float64}, params[:order])
    reffeₚ = ReferenceFE(lagrangian, Float64, params[:order])


    V = TestFESpace(params[:model], reffeᵤ, conformity=:H1, dirichlet_tags=u_diri_tags)
    U = TransientTrialFESpace(V, u_diri_values)

    Q = TestFESpace(params[:model], reffeₚ, conformity=:H1, dirichlet_tags=p_diri_tags)
    P = TrialFESpace(Q, p_diri_values)

    Y = MultiFieldFESpace([V, Q])
    X = TransientMultiFieldFESpace([U, P])

    return V, U, P, Q, Y, X
end
"""
  create_initial_conditions(params::Dict{Symbol,Any})

It creates the initial conditions for velocity and pressure. If `restart` is `true` then the velocity and the pressure field are interpoled on the specified DataFrame.  
"""
function create_initial_conditions(params::Dict{Symbol,Any})
    @unpack U,P,u0, D, restart,t0,Ω, printinitial = params


        uh0 = interpolate_everywhere(u0(t0), U(t0))
        ph0 = interpolate_everywhere(0.0, P(t0))

        if !restart
          if haskey(params,:p0)
            @unpack p0 = params
            ph0 = interpolate_everywhere(p0(t0), P(t0))
          
          end
        else

          tree = create_search_tree(params)
          uh_0 = restart_uh_field(params,tree)
          ph_0 = restart_ph_field(params,tree)
          uh0 = interpolate_everywhere(uh_0, U(t0))
          ph0 = interpolate_everywhere(ph_0, P(t0))

        end

        if printinitial
          writevtk(Ω,"Initial_Conditions_", cellfields=["uh0"=>uh0,"ph0"=>ph0])
        end
  
  return uh0,ph0
end

"""
  create_PETSc_setup(M::AbstractMatrix,ksp_setup::Function)

Wrapper for creating PETSc symbolic and numeric setup for `GridapPETSc`  
"""
function create_PETSc_setup(M::AbstractMatrix,ksp_setup::Function)
      solver = PETScLinearSolver(ksp_setup)
      ss = symbolic_setup(solver, M)
      ns = numerical_setup(ss, M)
      @check_error_code GridapPETSc.PETSC.KSPView(ns.ksp[],C_NULL)

      return ns
end

"""
  solve_case(params::Dict{Symbol,Any})

It solves iteratively the velocity and pressure system.
"""
function solve_case(params::Dict{Symbol,Any})

create_export_tags!(params)

@unpack M, petsc_options, time_step, θ, dt,t0, case, benchmark, method, trials, tests, 
Ω, matrix_freq_update, a_err_threshold, log_dir, save_sim_dir = params
@unpack U,P,u0 = params


uh0, ph0 = create_initial_conditions(params)

matrices = initialize_matrices_and_vectors(trials,tests, t0+dt, uh0, params; method=method)

Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, 
Mat_App, Mat_ML, Mat_inv_ML, Mat_S, Vec_Au, Vec_Ap = matrices



vec_pm,vec_um,vec_am,vec_sum_pm,Δa_star,Δpm1,Δa,b1,b2,ũ_vector = initialize_vectors(matrices,uh0,ph0)

if case == "TaylorGreen"
  @unpack u0,p0 = params
end

get_local_unique_idx(params)
export_nodes_glob(params)
export_n_Γ(params)

mkpath(save_sim_dir)



uh_avg = FEFunction(U(t0), vec_um)
set_zeros!(uh_avg.fields)
# uh_avg.fields.item.free_values .= 0.0
ph_avg = FEFunction(P(t0), vec_pm)
# ph_avg.fields.item.free_values .= 0.0
set_zeros!(ph_avg.fields)

GridapPETSc.with(args=split(petsc_options)) do
ns1 = create_PETSc_setup(Mat_ML,vel_kspsetup)
ns2 = create_PETSc_setup(Mat_S,pres_kspsetup)
uh_tn_updt = FEFunction(U(t0), vec_um)

for (ntime,tn) in enumerate(time_step)

      m = 0

      if mod(ntime,matrix_freq_update)==0
        println("update_matrices")
    
        @time begin
         Mat_Tuu, Mat_Tpu, Mat_Auu, Mat_Aup, Mat_Apu, Mat_App, 
          Mat_ML, Mat_inv_ML, Mat_S, Vec_Au, Vec_Ap = matrices_and_vectors(trials, tests, tn+dt, uh_tn_updt, params; method=method)
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
        
        M
      
      while (m<= M) && (err_norm_Δa0<a_err_threshold)

        Δpm1 .=  pazeros(Mat_S)
        Δa_star .= pazeros(Mat_ML)

        vec_um = GridapDistributed.change_ghost(vec_um, Mat_Auu)
        vec_pm = GridapDistributed.change_ghost(vec_pm, Mat_Aup)
        vec_am = GridapDistributed.change_ghost(vec_am, Mat_ML)


        println("solving velocity")
          
          b1 .= -Mat_Auu * vec_um - Mat_Aup * vec_pm - Mat_ML * vec_am +
          Mat_Auu * dt * vec_am + (1 - θ) * Mat_Aup * vec_sum_pm + Vec_Au

          @time solve!(Δa_star,ns1,b1)
          


        vec_um = GridapDistributed.change_ghost(vec_um, Mat_Apu)
        vec_pm = GridapDistributed.change_ghost(vec_pm, Mat_App)
        vec_am = GridapDistributed.change_ghost(vec_am, Mat_Tpu)

        Δa_star = GridapDistributed.change_ghost(Δa_star, Mat_Tpu)

        println("solving pressure")

        #-Vec_A because changing sign in the continuity equations
        b2 .= Mat_Tpu * Δa_star + Mat_Apu * (vec_um + dt * Δa_star) + Mat_App * vec_pm + Mat_Tpu * vec_am - Vec_Ap

        @time solve!(Δpm1,ns2,b2)

      
      println("update end")
      Δpm1 = GridapDistributed.change_ghost(Δpm1, Mat_Aup)

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

        println("err a")
        println(err_norm_Δa0)

        println("err p")
        println(err_norm_Δp0)

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
save_path = joinpath(save_sim_dir,"$(case)_$(tn)_.vtu")
update_time_average!(uh_tn,ph_tn, uh_avg, ph_avg, tn, params)

  if !benchmark && (mod(ntime,100)==0 || print_on_request(log_dir)) 


    if case == "TaylorGreen"
        writevtk(Ω, save_path, cellfields = ["uh" => uh_tn, "uh_analytic"=> u0(tn), "ph" => ph_tn, "ph_analytic"=> p0(tn)])
    else
      @time writevtk(Ω, save_path, cellfields = ["uh" => uh_tn, "uh_updt" => uh_tn_updt, "ph" => ph_tn,
      "uh_avg" => uh_avg,"ph_avg" => ph_avg])

    end
  end
  
  # read_nodes_to_export(params,uh_tn,tn)
  export_fields(params, tn, uh_tn, ph_tn)



  end #end for
end #end GridapPETSc


end #end solve_case