"""
  cconv(uadv, ∇u) 

Wrapper for the convective term
  ``u(\\nabla u)``  
"""
cconv(u, ∇u) = (∇u') ⋅ u

function segregated_equations(u_adv,params::Dict{Symbol,Any},simcase::SimulationCase)
  @sunpack skew, ν,dt, θ,D, order = simcase
  
  sprob = simcase.sprob
  @unpack dΩ, Ω = params
  @unpack skew = sprob
    
    skewcoeff = skew * 0.5 # ==0 if skew == false
    vms_activation = is_VMS(sprob.method)
    
    stab_coeff = compute_stab_coeff(simcase,params)
    Tm = momentum_stabilization(u_adv, stab_coeff, simcase)
    Tc = continuity_stabilization(u_adv, stab_coeff, simcase)
  
    Tuu(u, v) = ∫((v + (Tm*u_adv) ⋅ (∇(v)+(∇(v))')) ⊙ u)dΩ
    Tpu(u, q) = ∫(Tm * (∇(q)) ⊙ u)dΩ

    function Auu(u, v)
      uconv = u_adv ⋅ ∇(u)
      vconv = u_adv ⋅ ∇(v)
      return ∫(ν * ∇(v) ⊙ ∇(u) + uconv ⋅ v + (Tm * (vconv + u_adv ⋅ (∇(v))')) ⊙ uconv + (Tc * (∇ ⋅ v)) ⊙ (∇ ⋅ u))dΩ
    end
  
    Aup(p, v) = ∫((v + (Tm*u_adv) ⋅ (∇(v)+(∇(v))')) ⊙ ∇(p))dΩ
  
    Apu(u, q) = ∫(q * (∇ ⋅ u) + Tm * (∇(q)) ⊙ (u_adv ⋅ ∇(u)))dΩ
  
    App(p, q) = ∫((Tm * ∇(q)) ⊙ (∇(p)))dΩ
  
    S(p, q) = - θ * ∫((dt .+ Tm) ⋅ ((∇(q))' ⊙ (∇(p))))dΩ

    rhs(v) = 0.0

    return Tuu,Tpu,Auu,Aup,Apu,App,S,rhs

end

function is_VMS(method::VMS)
  return 1
end

function is_VMS(method::SUPG)
  return 0
end

    

