using PartitionedArrays
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using MPI
using Test


function TGV_Natural_test(backend)
    TGV_Natural_test(backend,2)
    TGV_Natural_test(backend,3)
end


function TGV_Natural_test(backend, D::Int64)

    t0 =0.0
    dt = 1e-2
    tF = dt * 3
    vortex_diameter = 1.0
    
    N = 32
    
    Re = 500
    rank_partition =  ntuple(i -> 2, D)


    sprob = StabilizedProblem(VMS(1))
    timep = TimeParameters(t0=t0,dt=dt,tF=tF)

    physicalp = PhysicalParameters(Re=Re,c=vortex_diameter)
    solverp = SolverParameters(matrix_freq_update = 1, M=2)
    exportp = ExportParameters(printinitial=true,printmodel=true)


    meshp= MeshParameters(rank_partition,D;N=N,L=vortex_diameter/2)

    simparams = SimulationParameters(timep,physicalp,solverp,exportp)

    bc_tgv = Natural(meshp,physicalp) 


    mcase = TaylorGreen(bc_tgv, meshp,simparams,sprob)



    @test SegregatedVMSSolver.solve(mcase,backend)

end


