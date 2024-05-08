module StabilizedEquationsTests

using Test
using SegregatedVMSSolver
using Gridap
using GridapDistributed
using PartitionedArrays

using SegregatedVMSSolver.Equations



function test_stab(rank_partition, distribute, D)
    ν = 0.001
    dt = 0.1
    θ = 0.5
    order = 1

    parts = distribute(LinearIndices((prod(rank_partition),)))
    domain = (D == 2) ? (0, 1, 0, 1) : (0, 1, 0, 1, 0, 1)
    mesh_partition = (D == 2) ? (4, 4) : (4, 4, 4)
    model = CartesianDiscreteModel(parts, rank_partition, domain, mesh_partition)
    Ω = Triangulation(model)
    dΩ = Measure(Ω, 2 * order)

    
    sprob = StabilizedProblem()

    params = Dict(:D => D, :ν => ν, :dt => dt, :θ => θ, :dΩ => dΩ, :Ω => Ω, :sprob=>sprob)

    vz((x, y)) = (D == 2) ? VectorValue(0.0, 0.0) : VectorValue(0.0, 0.0, 0.0)
    va((x, y)) = (D == 2) ? VectorValue(1.0, 2.0) : VectorValue(1.0, 2.0, 3.0)

    reffe = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)
    V = TestFESpace(model, reffe, dirichlet_tags="boundary")
    U = TrialFESpace(vz, V)
    uadv = interpolate_everywhere(va, U)

    vms_equations = segregated_equations(uadv, params)

    #Test creation VMS equations
    @testset "VMS Segregated Equations" for f in vms_equations
        @test typeof(f) <: Function
    end

    sprob = StabilizedProblem(SUPG(), ScalarFormulation(), true)
    params[:sprob] = sprob

    #Test creation SUPG equations
    supg_equations = segregated_equations(uadv, params)
    @testset "SUPG Segregated Equations" for f in supg_equations
        @test typeof(f) <: Function
    end
end

function main(distribute)
    for D in [2, 3]
        rank_partition = (D == 2) ? (2, 2) : (2, 2, 1)
        test_stab(rank_partition, distribute, D)
    end
end

end


