
abstract type StabilizationMethod end
abstract type StabilizationFormulation end
abstract type StabilizationParameters end


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
    order::Int64
end

struct SUPG <: StabilizationMethod
    order::Int64
end

function StabilizedProblem()
    StabilizedProblem(VMS(1))
end

function StabilizedProblem(method::VMS)
    StabilizedProblem(method, TensorFormulation(), false)
end

function StabilizedProblem(method::VMS, vv::Vector{Int64})
    StabilizedProblem(method, TensorFormulation(vv), false)
end

function StabilizedProblem(method::SUPG)
    StabilizedProblem(method, ScalarFormulation(), true)
end


#Default constructurs
function TensorFormulation()
    TensorFormulation(2,[4,36])
end

function TensorFormulation(vv::Vector{Int64})
    TensorFormulation(2,vv)
end

function ScalarFormulation()
    ScalarFormulation(1)
end