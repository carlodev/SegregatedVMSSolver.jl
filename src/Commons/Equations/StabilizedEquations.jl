"""
  cconv(uadv, ∇u) 

Wrapper for the convective term
  ``u(\\nabla u)``  
"""
cconv(u, ∇u) = (∇u') ⋅ u

function segregated_equations(u_adv,params::Dict{Symbol,Any},simcase::SimulationCase)
  @sunpack skew, ν,dt, θ,D = simcase
  
  sprob = simcase.sprob
  @unpack dΩ = params
  @unpack skew = sprob
    
    skewcoeff = skew * 0.5 # ==0 if skew == false
    vms_activation = is_VMS(sprob.method)
    
    visc_term =   0.0 #- ν* compute_Δ(u_adv,D) # compute the laplacian of the previous time step
    stab_coeff = compute_stab_coeff(simcase,params)
    Tm = momentum_stabilization(u_adv, stab_coeff, simcase)
    Tc = continuity_stabilization(u_adv, stab_coeff, simcase)
  
    ### VMS extra equations
    Auu_vms1(u, v) =  ∫(vms_activation .* u_adv ⋅ (∇(v))'*Tm⊙ (cconv ∘ (u_adv, ∇(u)) ) )dΩ + 
     ∫((vms_activation *skewcoeff) .* u_adv ⋅ (∇(v))'*Tm⊙ (u_adv ⋅ (∇ ⋅ u)) )dΩ

    Tuu_vms1(u, v) =  ∫(vms_activation .* u_adv ⋅ (∇(v))'*Tm⊙u )dΩ # VMS1

    Aup_vms1(p, v) =  ∫(vms_activation .* u_adv ⋅ (∇(v))'*Tm⊙ (∇(p)) )dΩ 

    #### SUPG equations
    Tuu(u, v) = ∫((v + Tm * (cconv ∘ (u_adv, ∇(v)))) ⊙ u)dΩ + Tuu_vms1(u, v) 
    Tpu(u, q) = ∫(Tm * (∇(q)) ⊙ u)dΩ
    Auu1(u, v) = ∫(ν * ∇(v) ⊙ ∇(u) + (cconv ∘ (u_adv, ∇(u))) ⋅ v + (Tm * (cconv ∘ (u_adv, ∇(v)))) ⊙ (cconv ∘ (u_adv, ∇(u)) ))dΩ
  
    Auu2(u, v) = ∫((Tc * (∇ ⋅ v)) ⊙ (∇ ⋅ u) + skewcoeff .* u_adv ⋅ (v + Tm * (cconv ∘ (u_adv, ∇(v)))) ⋅ (∇ ⋅ u))dΩ
  
    Auu(u, v) = Auu1(u, v) + Auu2(u, v) + Auu_vms1(u, v)
  
    Aup(p, v) = ∫( v⋅(∇(p) ) + (Tm * (cconv ∘ (u_adv, ∇(v)))) ⊙ ∇(p))dΩ + Aup_vms1(p, v)
  
    Apu(u, q) = ∫(q * (∇ ⋅ u) + skewcoeff .* Tm ⋅ (∇(q)) ⋅ u_adv ⋅ (∇ ⋅ u) + Tm * (∇(q)) ⊙ (cconv ∘ (u_adv, ∇(u)) ))dΩ
  
    App(p, q) = ∫((Tm * ∇(q)) ⊙ (∇(p)))dΩ
  
    ML(u, v) = Tuu(u, v) + (θ * dt) * Auu(u, v)
  
    S(p, q) = - θ * ∫((dt .+ Tm) ⋅ ((∇(q))' ⊙ (∇(p))))dΩ

    rhs(v) = 0.0
    
    return Tuu,Tpu,Auu,Aup,Apu,App,ML,S,rhs

end

function is_VMS(method::VMS)
  return 1
end

function is_VMS(method::SUPG)
  return 0
end


