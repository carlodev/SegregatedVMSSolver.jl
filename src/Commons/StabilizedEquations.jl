
cconv(uadv, ∇u) = uadv ⋅ (∇u)

val(x) = x
val(x::Gridap.Fields.ForwardDiff.Dual) = x.value 




function segregated_equations_SUPG!(u_adv, params)
    @unpack ν, dt, dΩ, D, Ω, θ = params

    h = h_param(Ω, D)
    merge!(params, Dict(:h => h))

     function τsu(u, h)
      r = 1
      τ₂ = h^2 / (4 * ν)
      τ₃ = dt / 2
    
      u = val(norm(u))
      if iszero(u)
        return (1 / τ₂^r + 1 / τ₃^r)^(-1 / r)
      end
      τ₁ = h / (2 * u)
      return (1 / τ₁^r + 1 / τ₂^r + 1 / τ₃^r)^(-1 / r)
    
    end
    
    function τb(u, h)
      return (u ⋅ u) * τsu(u, h)
    end
    
    
    Tuu(u, v) = ∫((v + τsu ∘ (u_adv, h) * (cconv ∘ (u_adv, ∇(v)))) ⊙ u)dΩ
    Tpu(u, q) = ∫((τsu ∘ (u_adv, h)) * (∇(q)) ⊙ u)dΩ
    Auu1(u, v) = ∫(ν * ∇(v) ⊙ ∇(u) + (cconv ∘ (u_adv, ∇(u))) ⋅ v + ((τsu ∘ (u_adv, h)) * (cconv ∘ (u_adv, ∇(v)))) ⊙ (cconv ∘ (u_adv, ∇(u))))dΩ
  
    Auu2(u, v) = ∫(((τb ∘ (u_adv, h)) * (∇ ⋅ v)) ⊙ (∇ ⋅ u) + 0.5 .* u_adv ⋅ (v + (τsu ∘ (u_adv, h)) * (cconv ∘ (u_adv, ∇(v)))) ⋅ (∇ ⋅ u))dΩ
  
    Auu(u, v) = Auu1(u, v) + Auu2(u, v)
  
    Aup(p, v) = ∫(-(∇ ⋅ v) * p + ((τsu ∘ (u_adv, h)) * (cconv ∘ (u_adv, ∇(v)))) ⊙ ∇(p))dΩ
  
    Apu(u, q) = ∫(q * (∇ ⋅ u) + 0.5 .* (τsu ∘ (u_adv, h)) ⋅ (∇(q)) ⋅ u_adv ⋅ (∇ ⋅ u) + (τsu ∘ (u_adv, h)) * (∇(q)) ⊙ (cconv ∘ (u_adv, ∇(u))))dΩ
  
    App(p, q) = ∫(((τsu ∘ (u_adv, h)) * ∇(q)) ⊙ (∇(p)))dΩ
  
    ML(u, v) = Tuu(u, v) + (θ * dt) * Auu(u, v)
  
    S(p, q) = - θ * ∫((dt .+ τsu ∘ (u_adv, h)) ⋅ ((∇(q)') ⊙ (∇(p))))dΩ

    rhs(v) = 0.0

    return Tuu,Tpu,Auu,Aup,Apu,App,ML,S,rhs
  end
    


 function segregated_equations_VMS!(u_adv, params)
    @unpack ν, dt, dΩ, θ, Ω, Cᵢ = params

    G, GG, gg = G_params(Ω, params)
    merge!(params, Dict(:G => G, :GG => GG, :gg => gg))

    function τm(uu, G, GG)
      τ₁ = Cᵢ[1] * (2 / dt)^2 #Here, you can increse the 2 if CFL high
      τ₃ = Cᵢ[2] * (ν^2 * GG)
  
  
      D = length(uu)
      if D == 2
        uu1 = val(uu[1])
        uu2 = val(uu[2])
        uu_new = VectorValue(uu1, uu2)
      elseif D == 3
        uu1 = val(uu[1])
        uu2 = val(uu[2])
        uu3 = val(uu[3])
        uu_new = VectorValue(uu1, uu2, uu3)
      end
  
      if iszero(norm(uu_new))
        return (τ₁ .+ τ₃) .^ (-1 / 2)
      end
  
      τ₂ = uu_new ⋅ G ⋅ uu_new
      return (τ₁ .+ τ₂ .+ τ₃) .^ (-1 / 2)
    end
  
    function τc(uu, gg, G, GG)
      return 1 / (τm(uu, G, GG) ⋅ gg)
    end


    Tm = τm∘(u_adv, G, GG) 
    Tc = τc∘(u_adv, gg, G, GG)

    Tuu(u, v) = ∫(v ⋅ u)dΩ + ∫(u_adv ⋅ ∇(v)*Tm⊙u )dΩ + ∫(u_adv ⋅ (∇(v))'*Tm⊙u )dΩ #Galerkin, SUPG, VMS1
    Tpu(u, q) = ∫(Tm * (∇(q)) ⊙ u)dΩ

    Auu_g(u, v) = ∫(ν * ∇(v) ⊙ ∇(u) + (cconv ∘ (u_adv, ∇(u))) ⋅ v )dΩ 
    Auu_supg(u, v) = ∫((Tm * (cconv ∘ (u_adv, ∇(v)))) ⊙ (cconv ∘ (u_adv, ∇(u))))dΩ +∫((Tc * (∇ ⋅ v)) ⊙ (∇ ⋅ u))dΩ
    Auu_vms1(u, v) = ∫(u_adv ⋅ (∇(v))'*Tm⊙ (cconv ∘ (u_adv, ∇(u))) )dΩ
  
    Auu(u, v) = Auu_g(u, v) + Auu_supg(u, v) + Auu_vms1(u, v)
  
   
    Aup(p, v) = ∫(-(∇ ⋅ v) * p + (Tm * (cconv ∘ (u_adv, ∇(v)))) ⊙ ∇(p))dΩ +∫(u_adv ⋅ (∇(v))'*Tm⊙ (∇(p)) )dΩ #Galerkin, SUPG, VMS1
  
    Apu(u, q) = ∫(q * (∇ ⋅ u) + Tm*∇(q)⋅(cconv ∘ (u_adv, ∇(u))))dΩ #Galerkin, SUPG
  
    App(p, q) = ∫((Tm * ∇(q)) ⊙ (∇(p)))dΩ
  
    ML(u, v) = Tuu(u, v) + (θ * dt) * Auu(u, v)
  
    S(p, q) =  - θ * ∫((dt .+ Tm) ⋅ ((∇(q)') ⊙ (∇(p))))dΩ

    rhs(v) = 0.0

    return Tuu,Tpu,Auu,Aup,Apu,App,ML,S,rhs
  
  end

    