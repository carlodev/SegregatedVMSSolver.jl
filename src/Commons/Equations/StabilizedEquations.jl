"""
  cconv(uadv, ∇u) 

Wrapper for the convective term
  ``u(\\nabla u)``  
"""
cconv(u, ∇u) = (∇u') ⋅ u


function Rm_adv_update(un,un_1,pn,ν,dt, U, method)
  _, vms_cross_activation = VMS_activation(method)
  if vms_cross_activation 
    Rm = (un-un_1) ./ dt +  cconv(un, ∇(un)) + ∇(pn) - ν*Δ(un)
    return interpolate(Rm,U)
  else
    return nothing
  end
end


macro create_equation(name, vms_flag, vms_cross_flag, base_term, vms_term, vms_cross_term)
  func_name = esc(name)

  quote
      function $(func_name)(a, b)
          if $(esc(vms_flag)) && $(esc(vms_cross_flag))
              $(esc(base_term))(a, b) + $(esc(vms_term))(a, b) + $(esc(vms_cross_term))(a, b)
          elseif $(esc(vms_flag))
              $(esc(base_term))(a, b) + $(esc(vms_term))(a, b)
          else
              $(esc(base_term))(a, b)
          end
      end
  end
end



function segregated_equations(u_adv,params::Dict{Symbol,Any},simcase::SimulationCase)
  @sunpack skew, ν,dt, θ,D = simcase
  
    sprob = simcase.sprob
    @unpack dΩ = params
    @unpack skew = sprob
      
    skewcoeff = skew * 0.5 # ==0 if skew == false
    vms_activation, vms_cross_activation = VMS_activation(sprob.method)
    
    stab_coeff = compute_stab_coeff(simcase,params)
    Tm = momentum_stabilization(u_adv, stab_coeff, simcase)
    Tc = continuity_stabilization(u_adv, stab_coeff, simcase)
    

      if vms_cross_activation
        @unpack Rm_adv = params
      end

    ### VMS EXTRA TERMS
    Auu_vms1(u, v) =  ∫(u_adv ⋅ (∇(v))'*Tm⊙ (cconv ∘ (u_adv, ∇(u)) ) )dΩ + 
    ∫( skewcoeff .* u_adv ⋅ (∇(v))'*Tm⊙ (u_adv ⋅ (∇ ⋅ u)) )dΩ
    Tuu_vms1(u, v) =  ∫( u_adv ⋅ (∇(v))'*Tm⊙u )dΩ # VMS1
    Aup_vms1(p, v) =  ∫( u_adv ⋅ (∇(v))'*Tm⊙ (∇(p)) )dΩ 

    ### VMS CROSS TERMS
    TmRm = Tm .*Rm_adv

    Tuu_vms_cross(u, v) = -0.5 * ∫( ((∇(v))⊙outer(TmRm, Tm⊙ u)) + ((∇(v)) ⊙outer(Tm⊙ u, TmRm)) )dΩ
    Auu_vms_cross(u, v) = -0.5* ∫((∇(v)⊙outer(TmRm, Tm⊙ (cconv ∘ (u_adv, ∇(u))))) -(∇(v)⊙outer(Tm⊙ (cconv ∘ (u_adv, ∇(u))), TmRm)) )dΩ
    Aup_vms_cross(p, v) = -0.5* ∫((∇(v)⊙outer(TmRm, Tm⊙ ∇(p))) +(∇(v)⊙outer(Tm⊙ ∇(p), TmRm)) )dΩ

    
    

    #### SUPG+GALERKIN TERMS
    Tuu_G_SUPG(u, v) = ∫((v + Tm * (cconv ∘ (u_adv, ∇(v)))) ⊙ u)dΩ + Tuu_vms1(u, v) + Tuu_vms_cross(u, v) 

    Auu_G_SUPG(u, v) = ∫(ν * ∇(v) ⊙ ∇(u) + (cconv ∘ (u_adv, ∇(u))) ⋅ v + (Tm * (cconv ∘ (u_adv, ∇(v)))) ⊙ (cconv ∘ (u_adv, ∇(u)) ))dΩ+ 
        ∫((Tc * (∇ ⋅ v)) ⊙ (∇ ⋅ u) + skewcoeff .* u_adv ⋅ (v + Tm * (cconv ∘ (u_adv, ∇(v)))) ⋅ (∇ ⋅ u))dΩ

    Aup_G_SUPG(p, v) = ∫( v⋅(∇(p) ) + (Tm * (cconv ∘ (u_adv, ∇(v)))) ⊙ ∇(p))dΩ + Aup_vms1(p, v) + Aup_vms_cross(p,v)

    @create_equation(Tuu, vms_activation, vms_cross_activation, Tuu_G_SUPG, Tuu_vms1, Tuu_vms_cross)
    @create_equation(Auu, vms_activation, vms_cross_activation, Auu_G_SUPG, Auu_vms1, Auu_vms_cross)
    @create_equation(Aup, vms_activation, vms_cross_activation, Aup_G_SUPG, Aup_vms1, Aup_vms_cross)



    Tpu(u, q) = ∫(Tm * (∇(q)) ⊙ u)dΩ

    Apu(u, q) = ∫(q * (∇ ⋅ u) + skewcoeff .* Tm ⋅ (∇(q)) ⋅ u_adv ⋅ (∇ ⋅ u) + Tm * (∇(q)) ⊙ (cconv ∘ (u_adv, ∇(u)) ))dΩ
  
    App(p, q) = ∫((Tm * ∇(q)) ⊙ (∇(p)))dΩ
     
    ML(u, v) = Tuu(u, v) + (θ * dt) * Auu(u, v)
  
    S(p, q) = - θ * ∫((dt .+ Tm) ⋅ ((∇(q))' ⊙ (∇(p))))dΩ

    rhs(v) = 0.0
    
    return Tuu,Tpu,Auu,Aup,Apu,App,ML,S,rhs

end

function VMS_activation(method::VMS)
    return true, method.cross_terms
end


function VMS_activation(method::SUPG)
    return false, false
end



