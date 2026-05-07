val_u(x) = x
val_u(x::Gridap.Fields.ForwardDiff.Dual) = x.value

# ── Effective Δt computation ───────────────────────────────────────────────────

"""
    compute_dt_eff(dt, limiter, uun, G, GG, ν) -> Float64

Return the effective Δt to be used in τ_M, based on the chosen limiter.
All methods guarantee Δt_eff ≥ dt.
"""
function compute_dt_eff(dt::Float64, ::NoLimiter, uun, G, GG, ν::Float64)
    return dt
end

function compute_dt_eff(dt::Float64, limiter::AdvectiveLimiter, uun, G, GG, ν::Float64)
    uu_new = VectorValue(val_u.(uun)...)
    iszero(norm(uu_new)) && return dt
    adv = uu_new ⋅ G ⋅ uu_new
    iszero(adv) && return dt
    return max(dt, limiter.C_adv / sqrt(adv))
end

function compute_dt_eff(dt::Float64, limiter::DiffusiveLimiter, uun, G, GG, ν::Float64)
    iszero(ν) && return dt
    return max(dt, limiter.C_diff / (ν * sqrt(GG)))
end

function compute_dt_eff(dt::Float64, limiter::CombinedLimiter, uun, G, GG, ν::Float64)
    uu_new  = VectorValue(val_u.(uun)...)
    adv     = uu_new ⋅ G ⋅ uu_new
    dt_adv  = (iszero(norm(uu_new)) || iszero(adv)) ? Inf : limiter.C_adv / sqrt(adv)
    dt_diff = iszero(ν) ? Inf : limiter.C_diff / (ν * sqrt(GG))
    return max(dt, min(dt_adv, dt_diff))
end

function compute_dt_eff(dt::Float64, limiter::CustomLimiter, uun, G, GG, ν::Float64)
    return limiter.f(dt, uun, G, GG, ν)
end

# ── Stabilization coefficient dispatch ────────────────────────────────────────

function compute_stab_coeff(coeff_method::ScalarFormulation, Ω, D::Int64)
    h = h_param(Ω, D)
    return ScalarStabilization(h)
end

function compute_stab_coeff(coeff_method::TensorFormulation, Ω, D::Int64)
    G, GG, gg = G_params(Ω, D)
    return TensorStabilization(G, GG, gg)
end

function compute_stab_coeff(simcase::SimulationCase, params::Dict{Symbol,Any})
    @unpack Ω = params
    @sunpack D, coeff_method = simcase
    compute_stab_coeff(coeff_method, Ω, D)
end

# ── Momentum stabilization ─────────────────────────────────────────────────────

"""
    momentum_stabilization(uu, stab_coeff::TensorStabilization, simcase)

τ_M from Bazilevs et al. (2007), with an optional Δt floor enforced
through the `dt_limiter` field of `TensorFormulation`. This prevents
τ_M → 0 (and consequently τ_C → ∞) when Δt is very small relative to
the mesh advective or diffusive time scales.
"""
function momentum_stabilization(uu, stab_coeff::TensorStabilization, simcase::SimulationCase)
    @unpack G, GG, gg = stab_coeff
    @unpack sprob = simcase
    @unpack Ci, dt_limiter = sprob.coeff_method

    @sunpack ν, dt = simcase

    function τm(uun, G, GG)
        dt_eff = compute_dt_eff(dt, dt_limiter, uun, G, GG, ν)
        println("Stabilization Using: dt = $(dt_eff), while physical dt = $dt")

        τ₁ = Ci[1] * (2 / dt_eff)^2
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
    continuity_stabilization(uu, stab_coeff::TensorStabilization, simcase)

τ_C = 1 / (τ_M · g·g). Inherits the Δt floor through τ_M automatically.
"""
function continuity_stabilization(uu, stab_coeff::TensorStabilization, simcase::SimulationCase)
    @unpack gg = stab_coeff
    return 1 / (momentum_stabilization(uu, stab_coeff, simcase) ⋅ gg)
end

# Scalar formulation unchanged
function momentum_stabilization(uu, stab_coeff::ScalarStabilization, simcase::SimulationCase)
    @unpack h = stab_coeff
    @unpack sprob = simcase
    @unpack r = sprob.coeff_method

    @sunpack ν, dt = simcase

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

function continuity_stabilization(uu, stab_coeff::ScalarStabilization, simcase::SimulationCase)
    return (uu ⋅ uu) * momentum_stabilization(uu, stab_coeff, simcase)
end