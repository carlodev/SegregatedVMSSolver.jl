module Interfaces

using SyntheticEddyMethod
using PartitionedArrays
    
export @sunpack
export evaluate_convergence
export pvnorm



"""
    VirtualBox(ylims::Tuple{Real,Real}, zlims::Tuple{Real,Real} ; σ=0.1)

Utility to create VirtualBox, specifing the ylims of the inlet, the zlims for the virtual plane, and the σ eddy dimension. Refer to the original documentation of the package SynteticEddyMethod for more info, and specifing more advanced settings
"""
function SyntheticEddyMethod.VirtualBox(ylims::Tuple{Real,Real}, zlims::Tuple{Real,Real} ; σ=0.1)
    
    ymin, ymax = ylims
    zmin, zmax = zlims
    
    @assert ymin<ymax
    @assert zmin<zmax

    x = collect(-15*σ:σ:15*σ)
    y = collect(ymin:σ:ymax)
    z = collect(zmin:σ:zmax)
    Vboxinfo = VirtualBox(x,y,z,σ)
    return Vboxinfo
end

function evaluate_convergence(ratio::Float64, field::String)
    if ratio<1
        @warn "Simulation not converging, $field is diverging, error =$ratio"
    else
        println("$field error = $ratio")
    end
end


function pvnorm(v)
    sqrt(sum(map(vi -> sum(abs2, vi), own_values(v))))
end


end