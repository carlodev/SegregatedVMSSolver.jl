module TurbParametersTests

using Test
using Revise
using SegregatedVMSSolver
using Parameters
using PartitionedArrays
using SyntheticEddyMethod

using SegregatedVMSSolver.ParametersDef

function main(distribute)
    TI = 0.001    
    D = 3
    meshfile = joinpath(@__DIR__, "..","..", "models", "sd7003s_3D_simple.msh")
    rank_partition=(2,2,1)

    #VirtualBox creation
    
    Vboxinfo = VirtualBox((-1.0,1.0), (0.0,0.2); σ=0.0125)
    physicalp = PhysicalParameters(Re=1000)
    Vboxinfo.N
    
    turbulencep= TurbulenceParameters(TI, Vboxinfo, physicalp)
    sprob = StabilizedProblem()

    timep = TimeParameters(0.0,0.1, 1.0)

    solverp = SolverParameters(M=2)
    exportp = ExportParameters(printmodel=false)
    meshp= MeshParameters(rank_partition,2,meshfile)

    simparams = SimulationParameters(timep,physicalp,turbulencep,solverp,exportp)
    tcase = Airfoil(meshp,simparams,sprob)

    @test typeof(tcase) <: SimulationCase

    eval_point = [ rand(3) for i = 1:1:   1000    ]
   
    
    tt = 0.1  
    @unpack TurbulenceInlet,Eddies,   Vboxinfo, Re_stress = tcase.simulationp.turbulencep
    Vboxinfo.σ
    
    for ep in eval_point
        compute_fluctuation(ep, tt, (TurbulenceInlet,Eddies, 1.0, Vboxinfo, Re_stress,D,1.0))
    end
    
    

    
    
    
end


end #end module