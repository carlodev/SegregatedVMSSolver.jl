using PartitionedArrays
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using MPI
using Gridap



    N = 128
    dt = 0.025

    order = 2
    t0 =0.0
    vortex_diameter = 1.0
    L = 0.5 * vortex_diameter
    
    u_conv = 0.2
    Re = 1_000

    D = 2

    CFL_eff = u_conv*N * dt * order / (2*L)

    backend = with_debug
    
    rank_partition = (4,4)
    tF = 0.1 #dt*10

      
    println("------------------------")
    println("CFL = $(CFL_eff)")
    println("------------------------")

    function petsc_options_cstm()
        return " -vel_ksp_type gmres -vel_pc_type gamg  -vel_ksp_rtol 1.e-8 -vel_ksp_converged_reason \
        -pres_ksp_type cg -pres_pc_type gamg  -pres_ksp_rtol 1.e-8 -pres_ksp_converged_reason"
    end


    my_lim = CustomLimiter() do dt, uun, G, GG, ν
        h = L /(order*N)
        dt_eff = h/u_conv
        return max(dt_eff,dt)

    end



    sprob = StabilizedProblem(VMS(order), TensorFormulation(dt_limiter = my_lim), false)




    # sprob = StabilizedProblem(method=VMS(order), coeff_method=TensorFormulation(r=2, Ci=[4, 36]), skew=false)
    timep = TimeParameters(t0=t0,dt=dt,tF=tF)

    physicalp = PhysicalParameters(Re=Re,c=vortex_diameter, u_in_mag=u_conv)

    solverp = SolverParameters(matrix_freq_update = 1, Number_Skip_Expansion=10e6, M = 40, a_err_threshold=10_000,
    petsc_options = petsc_options_cstm())
    exportp = ExportParameters(printinitial=true,printmodel=true, extra_export=["VelocityError","PressureError"])


    meshp= MeshParameters(rank_partition,D;N=N,L=L)


    simparams = SimulationParameters(timep,physicalp,solverp,exportp)

    params_tvg = TaylorGreenParameters(Vs=u_conv, Ua=0.3*u_conv, Va = 0.2*u_conv)
    bc_tgv = Periodic(meshp,physicalp,params_tvg ) 



    mcase = TaylorGreen(bc_tgv, meshp,simparams,sprob)


    SegregatedVMSSolver.solve(mcase,backend)




