
abstract type StabilizationMethod end
abstract type StabilizationFormulation end
abstract type StabilizationParameters end

"""
  cconv(uadv, ∇u) 

Wrapper for the convective term
  ``u(\\nabla u)``  
"""
conv(u, ∇u) = (∇u') ⋅ u


val(x) = x
val(x::Gridap.Fields.ForwardDiff.Dual) = x.value


struct StabilizedProblem{T<:StabilizationMethod,S<:StabilizationFormulation}
    method::T
    coeff_method::S
    skew::Bool
end



struct ScalarStabilization <: StabilizationParameters
    h
end

struct TensorStabilization <: StabilizationParameters
    G
    GG
    gg
end

struct ScalarFormulation <: StabilizationFormulation
    r::Int64
end

struct TensorFormulation <: StabilizationFormulation
    r::Int64
    Ci::Vector{Int64}
end


struct VMS <: StabilizationMethod
end

struct SUPG <: StabilizationMethod
end



function StabilizedProblem(method::VMS)
    StabilizedProblem(method, TensorFormulation(1,[4,36]), false)
end

function StabilizedProblem(method::VMS, vv::Vector{Int64})
    StabilizedProblem(method, TensorFormulation(1,vv), false)
end

function StabilizedProblem(method::SUPG)
    StabilizedProblem(method, ScalarFormulation(2), true)
end





function compute_stab_coeff(coeff_method::ScalarFormulation,params::Dict{Symbol,Any})
    @unpack Ω, D = params
    h = h_param(Ω, D,coeff_method)
    return ScalarStabilization(h)
end

function compute_stab_coeff(coeff_method::TensorFormulation,params::Dict{Symbol,Any})
    @unpack Ω = params
    G, GG, gg = G_params(Ω, params,coeff_method)
    return TensorStabilization(G, GG, gg)
end


function compute_stab_coeff(params::Dict{Symbol,Any})
    @unpack prob = params
    return compute_stab_coeff(prob.coeff_method, params)
end

function compute_stab_coeff(prob::StabilizedProblem,params::Dict{Symbol,Any})
    return compute_stab_coeff(prob.coeff_method, params)
end


"""
    momentum_stabilization(uu, stab_coeff::TensorStabilization, params::Dict{Symbol,Any} )

Stabilization parameter momentum stabilization
Bazilevs, Y., Calo, V. M., Cottrell, J. A., Hughes, T. J. R., Reali, A., & Scovazzi, G. (2007). Variational multiscale residual-based turbulence modeling for large eddy simulation of incompressible flows. Computer Methods in Applied Mechanics and Engineering, 197(1–4), 173–201. https://doi.org/10.1016/j.cma.2007.07.016
"""
function momentum_stabilization(uu, stab_coeff::TensorStabilization, params::Dict{Symbol,Any} )
    @unpack G,GG,gg=stab_coeff
    @unpack prob,ν,dt = params
    @unpack Ci = prob.coeff_method

    function τm(uu, G, GG)
        τ₁ = Ci[1] * (2 / dt)^2 #Here, you can increse the 2 if CFL high
        τ₃ = Ci[2] * (ν^2 * GG)




        uu_new = VectorValue(val.(uu)...)

        if iszero(norm(uu_new))
            return (τ₁ .+ τ₃) .^ (-1 / 2)
        end

        τ₂ = uu_new ⋅ G ⋅ uu_new
        return (τ₁ .+ τ₂ .+ τ₃) .^ (-1 / 2)
    end

    return τm ∘ (uu, G, GG)


end


"""
    continuity_stabilization(uu, stab_coeff::TensorStabilization, params::Dict{Symbol,Any} )

Stabilization parameter continuity
Bazilevs, Y., Calo, V. M., Cottrell, J. A., Hughes, T. J. R., Reali, A., & Scovazzi, G. (2007). Variational multiscale residual-based turbulence modeling for large eddy simulation of incompressible flows. Computer Methods in Applied Mechanics and Engineering, 197(1–4), 173–201. https://doi.org/10.1016/j.cma.2007.07.016
"""
function continuity_stabilization(uu, stab_coeff::TensorStabilization, params::Dict{Symbol,Any} )
     @unpack   gg = params
    return 1 / (momentum_stabilization(uu,stab_coeff,params) ⋅ gg)


end

"""
    momentum_stabilization(uu, stab_coeff::ScalarStabilization, params::Dict{Symbol,Any} )

Stabilization parameters momentum equation
Janssens, B. (2014). Numerical modeling and experimental investigation of ﬁne particle coagulation and dispersion in dilute ﬂows.
"""
function momentum_stabilization(uu, stab_coeff::ScalarStabilization, params::Dict{Symbol,Any} )
    @unpack h=stab_coeff
    @unpack prob,ν,dt = params
    @unpack r = prob.coeff_method

    function τsu(u, h)
        τ₂ = h^2 / (4 * ν)
        τ₃ = dt / 2

        u = val(norm(u))
        if iszero(u)
            return (1 / τ₂^r + 1 / τ₃^r)^(-1 / r)
        end
        τ₁ = h / (2 * u)
        return (1 / τ₁^r + 1 / τ₂^r + 1 / τ₃^r)^(-1 / r)

    end

    return τsu ∘ (uu, h)
end


"""
    continuity_stabilization(uu, stab_coeff::ScalarStabilization, params::Dict{Symbol,Any} )

Stabilization parameters continuity equation
Janssens, B. (2014). Numerical modeling and experimental investigation of ﬁne particle coagulation and dispersion in dilute ﬂows.
"""
function continuity_stabilization(uu, stab_coeff::ScalarStabilization, params::Dict{Symbol,Any} )
    @unpack   gg = params
    return (uu ⋅ uu) * momentum_stabilization(uu, stab_coeff, params)
end