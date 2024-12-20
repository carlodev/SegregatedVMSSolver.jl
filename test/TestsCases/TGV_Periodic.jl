using PartitionedArrays
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using MPI
using Test


function TGV_Periodic_test(backend)
    @test TGV_Periodic_test(backend,2)
    @test TGV_Periodic_test(backend,3)
end

function TGV_Periodic_test(backend,D::Int64)



    t0 =0.0
    dt = 0.01
    tF = dt * 3
    vortex_diameter = 1.0
    N = 32
    Re = 1600

    
    rank_partition = ntuple(i -> 2, D)


    


    solver_options = petsc_options(; vel_ksp="gmres", vel_pc="gamg", pres_ksp = "cg", pres_pc = "ilu")

    sprob = StabilizedProblem(VMS(1))
    timep = TimeParameters(t0=t0,dt=dt,tF=tF)

    physicalp = PhysicalParameters(Re=Re,c=vortex_diameter)
    solverp = SolverParameters(matrix_freq_update = 1, Number_Skip_Expansion=1, M = 40,
    petsc_options = solver_options)
    exportp = ExportParameters(printinitial=true,printmodel=true)


    meshp= MeshParameters(rank_partition,D;N=N,L= pi*vortex_diameter)


    simparams = SimulationParameters(timep,physicalp,solverp,exportp)

    bc_tgv = Periodic(meshp,physicalp ) 



    mcase = TaylorGreen(bc_tgv, meshp,simparams,sprob)


    @test SegregatedVMSSolver.solve(mcase,backend)

end
