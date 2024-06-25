module TurbParametersTests

using Test
using Revise
using SegregatedVMSSolver
using Parameters
using PartitionedArrays
using SyntheticEddyMethod

using SegregatedVMSSolver.ParametersDef

# function main(distribute)
    TI = 0.001    
    D = 3
    meshfile = joinpath(@__DIR__, "..", "models", "sd7003s_3D_simple.msh")
    rank_partition=(2,2,1)

    #VirtualBox creation
    Vboxinfo = VirtualBox((-6,6), (-0.50,0.50); σ=0.02)
    physicalp = PhysicalParameters(Re=1000)

    turbulencep= TurbulenceParameters(TI, Vboxinfo, physicalp)
    sprob = StabilizedProblem()

    timep = TimeParameters(0.0,0.1, 1.0)

    solverp = SolverParameters(M=2)
    exportp = ExportParameters(printmodel=false)
    meshp= MeshParameters(rank_partition,2,meshfile)

    simparams = SimulationParameters(timep,physicalp,turbulencep,solverp,exportp)
    tcase = Airfoil(meshp,simparams,sprob)

    @test typeof(tcase) <: SimulationCase

    eval_point = [0.0, 1.0, 0.0] #Point in the domain
    # for tt = 0.0:0.1:0.5
    tt = 0.1  
        u_fluct = compute_fluctuation(eval_point, tt, tcase)
        @test typeof(u_fluct) <: AbstractVector
    # end

end


end