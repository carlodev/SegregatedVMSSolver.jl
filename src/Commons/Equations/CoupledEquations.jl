function coupled_equations(params::Dict{Symbol,Any},simcase::SimulationCase)
    @sunpack skew, ν,dt, θ, D = simcase
    
    sprob = simcase.sprob
    @unpack dΩ = params
    @unpack skew = sprob

    #Conservations equations
    Rc(u) = ∇ ⋅ u
    dRc(du) = ∇ ⋅ du
    Rm(t, (u, p)) = ∂t(u) + u ⋅ ∇(u) + ∇(p)  #- ν*Δ(u)
    dRm((u, p), (du, dp), (v, q)) = du ⋅ ∇(u) + u ⋅ ∇(du) + ∇(dp) #- ν*Δ(du)
    stab_coeff = compute_stab_coeff(simcase,params)

    Tm(u) = momentum_stabilization(u, stab_coeff, simcase)
    Tc(u) = continuity_stabilization(u, stab_coeff, simcase)
  

    θvp = 1.0

    #VMS equations
    #Function used in VMS terms
    TRm(t, (u, p)) = Tm(u.cellfield) * Rm(t, (u, p))
  
    #Variational equations
    Bᴳ(t, (u, p), (v, q)) = ∫(∂t(u) ⋅ v)dΩ + ∫((u ⋅ ∇(u)) ⋅ v)dΩ - ∫((∇ ⋅ v) * p)dΩ + ∫((q * (∇ ⋅ u)))dΩ + ν * ∫(∇(v) ⊙ ∇(u))dΩ 
  
    #SUPG terms 
    B_SUPG(t, (u, p), (v, q)) = ∫((u ⋅ ∇(v) + ∇(q)) ⊙ TRm(t, (u, p)))dΩ + ∫((∇ ⋅ v) ⊙ (Tc(u.cellfield)  * Rc(u)))dΩ
  
    #First VMS term
    B_VMS1(t, (u, p), (v, q)) = ∫((u ⋅ (∇(v))') ⊙ TRm(t, (u, p)))dΩ
  
    #Second VMS term
    B_VMS2(t, (u, p), (v, q)) = -1 * ∫((∇(v)) ⊙ (outer(TRm(t, (u, p)), TRm(t, (u, p)))))dΩ
  
    #Adding all the contributions
    Bᴹ(t, (u, p), (v, q)) = Bᴳ(t, (u, p), (v, q)) + B_SUPG(t, (u, p), (v, q)) + B_VMS1(t, (u, p), (v, q)) + B_VMS2(t, (u, p), (v, q))
  
    res(t, (u, p), (v, q)) = Bᴹ(t, (u, p), (v, q))
  
    #VMS Jacobian
    #Function used in the VMS terms -  derivative
    dTRm((u, p), (du, dp), (v, q)) =  Tm(u.cellfield)* dRm((u, p), (du, dp), (v, q)) #Derivative of the Function used in the second VMS term
  
    #Variational equations - derivative
    dBᴳ(t, (u, p), (du, dp), (v, q)) = ∫(((du ⋅ ∇(u)) ⋅ v) + ((u ⋅ ∇(du)) ⋅ v) + (∇(dp) ⋅ v) + (q * (∇ ⋅ du)))dΩ + ν * ∫(∇(v) ⊙ ∇(du))dΩ
  
    #SUPG terms derivative
    dB_SUPG(t, (u, p), (du, dp), (v, q)) = ∫(((du ⋅ ∇(v)) ⋅ TRm(t, (u, p))) + ((u ⋅ ∇(v) + ∇(q)) ⋅ dTRm((u, p), (du, dp), (v, q))) + ((∇ ⋅ v) ⋅ ( Tc(u.cellfield)  .* dRc(du))))dΩ
  
    #First VMS term derivative
    dB_VMS1(t, (u, p), (du, dp), (v, q)) = ∫(((du ⋅ ∇(v)') ⊙ TRm(t, (u, p))) + ((u ⋅ ∇(v)') ⊙ dTRm((u, p), (du, dp), (v, q))))dΩ
  
  
    #Second VMS term derivative
    dB_VMS2(t, (u, p), (du, dp), (v, q)) = ∫(∇(v) ⊙ (outer(dTRm((u, p), (du, dp), (v, q)), TRm(t, (u, p)))))dΩ + ∫(∇(v) ⊙ (outer(TRm(t, (u, p)), dTRm((u, p), (du, dp), (v, q)))))dΩ
  
  
    #Adding all the contributions
    jac(t, (u, p), (du, dp), (v, q)) = dBᴳ(t, (u, p), (du, dp), (v, q)) + dB_SUPG(t, (u, p), (du, dp), (v, q)) + dB_VMS1(t, (u, p), (du, dp), (v, q)) - dB_VMS2(t, (u, p), (du, dp), (v, q))
  
    #VMS time-Jacobian
    dtTRm(t, (u, p), (dut, dpt), (v, q)) = ( Tm(u.cellfield)) ⋅ dut #Derivative of the Function used in the second VMS term
  
    dtBᴳ(t, (u, p), (dut, dpt), (v, q)) = ∫(dut ⋅ v)dΩ 
    dtB_SUPG(t, (u, p), (dut, dpt), (v, q)) =  ∫((u ⋅ ∇(v) + (θvp)*∇(q)) ⊙ dtTRm(t, (u, p), (dut, dpt), (v, q)))dΩ 
    dtB_VMS1(t, (u, p), (dut, dpt), (v, q)) = ∫((u ⋅ (∇(v))') ⊙ dtTRm(t, (u, p), (dut, dpt), (v, q)))dΩ
    dtB_VMS2(t, (u, p), (dut, dpt), (v, q)) =  ∫((∇(v)) ⊙ (outer(dtTRm(t, (u, p), (dut, dpt), (v, q)), TRm(t, (u, p)))))dΩ + ∫((∇(v)) ⊙ (outer(TRm(t, (u, p)), dtTRm(t, (u, p), (dut, dpt), (v, q)))))dΩ        
    
    jac_t(t, (u, p), (dut, dpt), (v, q)) = dtBᴳ(t, (u, p), (dut, dpt), (v, q)) +  dtB_SUPG(t, (u, p), (dut, dpt), (v, q)) + dtB_VMS1(t, (u, p), (dut, dpt), (v, q)) - dtB_VMS2(t, (u, p), (dut, dpt), (v, q))
  
    res, jac, jac_t
  
  end
  
      
  
  