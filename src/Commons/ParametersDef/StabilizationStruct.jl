
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
    τm_comp::Real = -1 #Standard VMS; τm_comp = 1 SUPG standard
    @assert length(Ci) == 2 "Ci length must be 2"
    @assert (Ci[1]>=0 && Ci[2]>=0) "Ci values must be non-negative"
    @assert τm_comp==1 ||  τm_comp==-1 " τm values +1 or -1 in TensorFormulation"

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

