
function intialize_picard_iteration(solver,Vg,V0,dΩ  )
    flux(∇u) = ∇u
    f(x) = 1
    a(u, v) = ∫(∇(v) ⊙ (flux ∘ ∇(u))) * dΩ
    l(v) = ∫( v * f) * dΩ

    op = AffineFEOperator(a, l, Vg, V0)

    uh = interpolate_everywhere(1.0, Vg)

    uh,cache = GridapPETSc.solve!(uh, solver, op)

    return uh
end


function run_picard_iteration(solver,Vg,V0,dΩ,uh_adv,p, γ)

    flux(∇u,∇uha) = norm(∇uha)^(p - 2) * ∇u
    f(x) = 1
    a(u, v) = ∫(∇(v) ⊙ (flux ∘ (∇(u),∇(uh_adv)))) * dΩ
    l(v) = ∫(v * f) * dΩ

    uh = interpolate_everywhere(1.0, Vg)

    op = AffineFEOperator(a, l, Vg, V0)
    uh,cache = GridapPETSc.solve!(uh, solver, op)
 

    return uh*γ + uh_adv*(1-γ)
end



function initial_velocity_wall_distance(params, walltag::String; pmax = 6)
    @unpack model, u_in, Re, c, D = params
    if pmax<6    @error "pmax is set $pmax, it has to be >5"  end

    options = "-ksp_type gmres -pc_type gamg  -ksp_monitor -snes_rtol 1.0e-10 -ksp_converged_reason -ksp_error_if_not_converged true"

  GridapPETSc.with(args=split(options)) do
      order = 1
      reffe = ReferenceFE(lagrangian, Float64, order)

      U0 = TestFESpace(model, ReferenceFE(lagrangian, VectorValue{D,Float64}, order); conformity=:H1, dirichlet_tags=walltag)
      zvector = (D==2) ? VectorValue(0.0, 0.0) : VectorValue(0.0,0.0, 0.0)
      Ug = TrialFESpace(U0, zvector)

      V0 = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags=walltag)
      Vg = TrialFESpace(V0, 0.0)

      degree = 2
      Ω = Triangulation(model)
      dΩ = Measure(Ω, degree)

      solver = PETScLinearSolver()

      uh2 = intialize_picard_iteration(solver,Vg,V0,dΩ)

      uh3 = run_picard_iteration(solver,Vg,V0,dΩ,uh2,3,1)
      uh = uh3
      
      for p in 4:1:pmax
        for inner_iter in 1:1:3
            uh = run_picard_iteration(solver,Vg,V0,dΩ,uh,p,0.5)
            println("p = $p, inner_iter = $(inner_iter)")
        end
      end
      


      p = pmax

    vnorm(up, ∇up) = -norm(∇up)^(p - 1) + (p / (p - 1) * up + norm(∇up)^(p))^((p - 1) / p)
    vh = vnorm ∘ (uh, ∇(uh))

    DistFeFun = interpolate_everywhere(vh, Vg)

   

    δ99 = bl_height(Re, c)

      function u_star_fun(d)
        xvel = bl_function(d,δ99,u_in)

    
          if D == 2
              return VectorValue(xvel, 0.0)
          else
              return VectorValue(xvel, 0.0, 0.0)

          end
      end


      u_start_oper = u_star_fun∘DistFeFun
      u_start_FEFun = interpolate_everywhere(u_start_oper,Ug)
      
      writevtk(Ω,"Inital_Condition", cellfields=["uh_start"=>u_start_FEFun, "vh"=>vh, "uh"=>uh, "uh3"=>uh3])
      
      u_start = get_free_dof_values(u_start_FEFun)
      
      return u_start

  end
end


function bl_height(Re::Real, L::Float64)
    if Re < 500e3
        δ99 = 4.91 * L / (Re^0.5)

    else
        δ99 = 0.38 * L / (Re^0.2)

    end
    return δ99

end


function bl_function(x,δ99::Float64,u_in::Float64)
    if x>δ99
      return u_in 
    else
    return (-(x/δ99) .^ 2 + 2 * (x/δ99))*u_in
    end
    
  end