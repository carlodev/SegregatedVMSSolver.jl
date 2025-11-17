module Projection
using Gridap
using GridapDistributed
using CSV, DelimitedFiles
using DataFrames
using Parameters
using PartitionedArrays
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.Interfaces
using SegregatedVMSSolver.CreateProblem


export compute_VMS2_error


###PROJECTION


function compute_VMS2_error(uh_fine, simcase::SimulationCase,params::Dict{Symbol,Any}, tn::Real)
    @unpack U, dΩ, degree,parts = params

    ubar, uprime = project_solution(uh_fine, simcase, params, tn)
    compute_stresses(ubar, uprime, dΩ) 

end

function create_coarse_spaces(params,simcase)
    compute_Uc = true
    haskey(params , Uc) ? compute_Uc = false : compute_Uc = true
    if compute_Uc
        @sunpack  model = params
        simcase_coarse = deepcapy(simcase)
        # simcase_coarse.meshp.meshinfo.N .= ones(Int64, D) .* Int(ceil(N[1] / 2))
        simcase_coarse.sprob.method.order = order -1
        boundary_conditions = create_boundary_conditions(simcase) 
        Vc, Uc, _, _ = creation_fe_spaces(simcase_coarse, model, boundary_conditions)
        merge!(params,Dict(:Uc=>Uc, :Vc=>Vc))
    end

    @unpack Vc, Uc = params

    return Vc, Uc
end

function project_solution(uh_fine, simcase::SimulationCase, params::Dict{Symbol,Any}, tn::Real)
return true
end

function project_solution(uh_fine, simcase::TaylorGreen{Periodic}, params::Dict{Symbol,Any}, tn::Real)
    @assert tn>=0.0
    @sunpack D,N, order = simcase
    @unpack U,V,P,Q, Ω, degree,parts, dΩ = params

    @assert order >1

    Vc, Uc =   create_coarse_spaces(params,simcase)

    #L2 projection of uh_fine solution on lower order dimensional space (same mesh)
    # Build weak form
    a(u,v) = ∫( u ⋅ v )dΩ
    l(v)   = ∫( uh_fine ⋅ v )dΩ

    # Assemble and solve
    op = AffineFEOperator(a,l,Uc(tn),Vc(tn))
    ubar = solve(op)

    println("ubar - solved")

    uprime = uh_fine - ubar

    return ubar, uprime
end

function compute_stresses(ubar, uprime, dΩ) # dΩ is of the coarse mesh
    tau_cross = ubar ⊗ uprime + uprime ⊗ ubar
    println("tau_cross")
    tau_re    = uprime ⊗ uprime
    println("tau_re")

    norm_cross = sum(∫( tau_cross⊙tau_cross )dΩ) #basically we compute the tensor product to compute the norm
    norm_re    = sum(∫( tau_re⊙tau_re )dΩ)

    R_stress = sqrt(norm_re ./ norm_cross)
    println("R_stress = $(R_stress)")

    grad_ubar = ∇(ubar)

    eps_cross = sum(∫( -(tau_cross ⊙ grad_ubar) )dΩ)
    eps_re    = sum(∫( -(tau_re⊙ grad_ubar) )dΩ)

        println("eps_cross = $(eps_cross)")
        println("eps_re = $(eps_re)")

end


end
