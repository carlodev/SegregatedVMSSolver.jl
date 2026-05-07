abstract type StabilizationMethod end
abstract type StabilizationFormulation end
abstract type StabilizationParameters end

# ── Timestep limiter hierarchy ─────────────────────────────────────────────────

abstract type TimestepLimiter end

"""
    NoLimiter()

Default behaviour: Δt is used as-is inside τ_M and τ_C.
No minimum is imposed on the temporal contribution 4/Δt².
"""
struct NoLimiter <: TimestepLimiter end

"""
    AdvectiveLimiter(; C_adv = 1.0)

Clamps Δt from below using the local advective time scale:

    Δt_eff = max(Δt, C_adv / √(u · G · u))

Prevents τ_M → 0 when Δt is small but the flow is fast relative to the mesh.
`C_adv` is an O(1) coefficient tuned so the limiter only activates well below
the advective CFL scale.
"""
@with_kw struct AdvectiveLimiter <: TimestepLimiter
    C_adv::Float64 = 1.0
    @assert C_adv > 0.0 "C_adv must be positive"
end

"""
    DiffusiveLimiter(; C_diff = 1.0)

Clamps Δt from below using the local diffusive time scale:

    Δt_eff = max(Δt, C_diff / (ν √(G:G)))

Useful for low-Reynolds or boundary-layer dominated flows where
diffusion sets the relevant time scale.
"""
@with_kw struct DiffusiveLimiter <: TimestepLimiter
    C_diff::Float64 = 1.0
    @assert C_diff > 0.0 "C_diff must be positive"
end

"""
    CombinedLimiter(; C_adv = 1.0, C_diff = 1.0)

Takes the smaller of the advective and diffusive time scales as the floor:

    Δt_adv  = C_adv  / √(u · G · u)
    Δt_diff = C_diff / (ν √(G:G))
    Δt_eff  = max(Δt, min(Δt_adv, Δt_diff))

Recommended for the Taylor-Green vortex validation and general LES
where both regimes may be encountered.
"""
@with_kw struct CombinedLimiter <: TimestepLimiter
    C_adv::Float64  = 1.0
    C_diff::Float64 = 1.0
    @assert C_adv  > 0.0 "C_adv must be positive"
    @assert C_diff > 0.0 "C_diff must be positive"
end

"""
    CustomLimiter(f)

User-supplied limiter. `f` must satisfy

    f(dt::Float64, uun, G, GG, ν::Float64) -> Float64

and return the effective Δt to use inside τ_M. The returned value must be
≥ `dt`; no safety check is performed, so the user is responsible for
ensuring stability.

# Example
```julia
my_lim = CustomLimiter((dt, uun, G, GG, ν) -> max(dt, 0.5 / (ν * sqrt(GG))))
```
"""
struct CustomLimiter <: TimestepLimiter
    f::Function
end

# ── Stabilization parameter containers ────────────────────────────────────────

@with_kw struct StabilizedProblem{T<:StabilizationMethod, S<:StabilizationFormulation}
    method::T           = VMS()
    coeff_method::S     = TensorFormulation()
    skew::Bool          = false
end

struct ScalarStabilization <: StabilizationParameters
    h
end

struct TensorStabilization <: StabilizationParameters
    G
    GG
    gg
end

@with_kw struct ScalarFormulation <: StabilizationFormulation
    r::Int64 = 1
end

@with_kw struct TensorFormulation <: StabilizationFormulation
    r::Int64                  = 2
    Ci::Vector{Real}          = [4, 36]
    dt_limiter::TimestepLimiter = NoLimiter()
    @assert (Ci[1] >= 0 && Ci[2] >= 0) "Ci values must be non-negative"
end

@with_kw struct VMS <: StabilizationMethod
    order::Int64 = 1
end

@with_kw struct SUPG <: StabilizationMethod
    order::Int64 = 1
end

function StabilizedProblem(method::VMS)
    StabilizedProblem(method, TensorFormulation(), false)
end

function StabilizedProblem(method::SUPG)
    StabilizedProblem(method, ScalarFormulation(), true)
end