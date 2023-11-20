module StabParamsTests
using Test
using SegregatedVMSSolver
using Gridap
using LinearAlgebra

N = 100
domains = [(0,1,0,1), (0,1,0,1,0,1)]
 for (i,D) in enumerate([2,3])

domain = domains[i]
partition = N .* ones(D)
model =  CartesianDiscreteModel(domain, partition)
Ω = Triangulation(model)
he = SegregatedVMSSolver.h_param(Ω, D)
@test norm(collect(he) .- 1/N .* ones(N^D)) <1e-10

G, GG, gg = SegregatedVMSSolver.G_params(Ω, Dict(:D=>D))
if D ==2
    @test norm( collect(G) .- TensorValue(N^D,0,0,N^D)) <1e-10
end #end if 
@test norm( collect(GG) .- D * N^4) <1e-10
@test norm( collect(gg) .- D .* N^2) <1e-10



end #end for 

end #end module 