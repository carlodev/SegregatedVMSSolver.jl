using PartitionedArrays
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using MPI
using Test


function TGV_Periodic_test()
    TGV_Periodic_test(2)
    TGV_Periodic_test(3)
end

function TGV_Periodic_test(D::Int64)

    t0 =0.0
    dt = 0.01
    tF = dt * 1
    vortex_diameter = 1.0
    N = 8
    Re = 1600

    
    rank_partition = ntuple(i -> 2, D)


    solver_options = petsc_options(; vel_ksp="gmres", vel_pc="gamg", pres_ksp = "cg", pres_pc = "gamg")

    sprob = StabilizedProblem(VMS(1))
    timep = TimeParameters(t0=t0,dt=dt,tF=tF)

    physicalp = PhysicalParameters(Re=Re,c=vortex_diameter)
    solverp = SolverParameters(matrix_freq_update = 1, Number_Skip_Expansion=10e10, M = 2,
    petsc_options = solver_options)
    exportp = ExportParameters(printinitial=false,printmodel=false)

    if D == 2
        l_mesh = 0.5
    elseif D== 3
        l_mesh = pi
    end

    meshp= MeshParameters(rank_partition,D;N=N,L= l_mesh*vortex_diameter)


    simparams = SimulationParameters(timep,physicalp,solverp,exportp)

    bc_tgv = Periodic(meshp,physicalp ) 



    mcase = TaylorGreen(bc_tgv, meshp,simparams,sprob)

    @test typeof(mcase) <: SimulationCase

    return mcase

end
