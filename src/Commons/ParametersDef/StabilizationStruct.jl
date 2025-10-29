
abstract type StabilizationMethod end
abstract type StabilizationFormulation end
abstract type StabilizationParameters end


@with_kw struct StabilizedProblem{T<:StabilizationMethod,S<:StabilizationFormulation}
    method::T = VMS()
    coeff_method::S=TensorFormulation()
    skew::Bool=false
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
    r::Int64 = 2 
    Ci::Vector{Real} = [4,36]
    @assert (Ci[1]>=0 && Ci[2]>=0) "Ci values must be non-negative"
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

