val_u(x) = x
val_u(x::Gridap.Fields.ForwardDiff.Dual) = x.value

function compute_stab_coeff(coeff_method::ScalarFormulation,Ω,D::Int64)
    h = h_param(Ω, D)
    return ScalarStabilization(h)
end

function compute_stab_coeff(coeff_method::TensorFormulation,Ω,D::Int64)
    G, GG, gg = G_params(Ω, D)
    return TensorStabilization(G, GG, gg)
end

function compute_stab_coeff(simcase::SimulationCase,params::Dict{Symbol,Any})
    @unpack Ω = params
    @sunpack D,coeff_method =simcase
    compute_stab_coeff(coeff_method,Ω,D)
end



"""
    momentum_stabilization(uu, stab_coeff::TensorStabilization,simcase::SimulationCase )

Stabilization parameter momentum stabilization
Bazilevs, Y., Calo, V. M., Cottrell, J. A., Hughes, T. J. R., Reali, A., & Scovazzi, G. (2007). Variational multiscale residual-based turbulence modeling for large eddy simulation of incompressible flows. Computer Methods in Applied Mechanics and Engineering, 197(1–4), 173–201. https://doi.org/10.1016/j.cma.2007.07.016
"""
function momentum_stabilization(uu, stab_coeff::TensorStabilization,simcase::SimulationCase )
    @unpack G,GG,gg = stab_coeff
    @unpack sprob = simcase
    @unpack Ci = sprob.coeff_method
    
    @sunpack ν, dt = simcase

    function τm(uun, G, GG)
        τ₁ = Ci[1] * (2 / dt)^2 #Here, you can increse the 2 if CFL high
        τ₃ = Ci[2] * (ν^2 * GG)

        uu_new = VectorValue(val_u.(uun)...)

        if iszero(norm(uu_new))
            return (τ₁ .+ τ₃) .^ (-1 / 2)
        end

        τ₂ = uu_new ⋅ G ⋅ uu_new
        return (τ₁ .+ τ₂ .+ τ₃) .^ (-1 / 2)
    end
 
    return τm ∘ (uu, G, GG)


end


"""
    continuity_stabilization(uu, stab_coeff::TensorStabilization,simcase::SimulationCase )

Stabilization parameter continuity
Bazilevs, Y., Calo, V. M., Cottrell, J. A., Hughes, T. J. R., Reali, A., & Scovazzi, G. (2007). Variational multiscale residual-based turbulence modeling for large eddy simulation of incompressible flows. Computer Methods in Applied Mechanics and Engineering, 197(1–4), 173–201. https://doi.org/10.1016/j.cma.2007.07.016
"""
function continuity_stabilization(uu, stab_coeff::TensorStabilization,simcase::SimulationCase)
     @unpack   gg = stab_coeff
    return 1 / (momentum_stabilization(uu,stab_coeff,simcase) ⋅ gg)


end

"""
    momentum_stabilization(uu, stab_coeff::ScalarStabilization,simcase::SimulationCase )

Stabilization parameters momentum equation
Janssens, B. (2014). Numerical modeling and experimental investigation of ﬁne particle coagulation and dispersion in dilute ﬂows.
"""
function momentum_stabilization(uu, stab_coeff::ScalarStabilization,simcase::SimulationCase )
    @unpack h=stab_coeff
    @unpack sprob = simcase
    @unpack r = sprob.coeff_method

    @sunpack ν,dt = simcase

    function τsu(u, h)
        τ₂ = h^2 / (4 * ν)
        τ₃ = dt / 2

        u = val_u(norm(u))
        if iszero(u)
            return (1 / τ₂^r + 1 / τ₃^r)^(-1 / r)
        end
        τ₁ = h / (2 * u)
        return (1 / τ₁^r + 1 / τ₂^r + 1 / τ₃^r)^(-1 / r)

    end

    return τsu ∘ (uu, h)
end


"""
    continuity_stabilization(uu, stab_coeff::ScalarStabilization,simcase::SimulationCase )

Stabilization parameters continuity equation
Janssens, B. (2014). Numerical modeling and experimental investigation of ﬁne particle coagulation and dispersion in dilute ﬂows.
"""
function continuity_stabilization(uu, stab_coeff::ScalarStabilization,simcase::SimulationCase )
    return (uu ⋅ uu) * momentum_stabilization(uu, stab_coeff, simcase)
end
