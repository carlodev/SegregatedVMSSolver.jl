using PartitionedArrays
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using MPI


function TGV_Periodic_test(backend)
    t0 =0.0
    dt = 0.02
    tF = dt * 5
    vortex_diameter = 1.0
    N = 16
    Re = 1000
    D = 2
    rank_partition = (2,2)


    sprob = StabilizedProblem(VMS(1))
    timep = TimeParameters(t0=t0,dt=dt,tF=tF)

    physicalp = PhysicalParameters(Re=Re,c=vortex_diameter)
    solverp = SolverParameters(matrix_freq_update = 1)
    exportp = ExportParameters(printinitial=false,printmodel=false)


    meshp= MeshParameters(rank_partition,D;N=16,L=vortex_diameter/2)


    simparams = SimulationParameters(timep,physicalp,solverp,exportp)

    bc_tgv = Periodic(meshp,physicalp ) 



    mcase = TaylorGreen(bc_tgv, meshp,simparams,sprob)


    @test SegregatedVMSSolver.solve(mcase,backend)

end
