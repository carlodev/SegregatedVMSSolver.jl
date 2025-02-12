using PartitionedArrays
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions

mesh_file = joinpath(@__DIR__,"..","..","models", "DU89_3D_A1_M.msh")



function petsc_options_cstm()
      return "-vel_ksp_type gmres -vel_ksp_gmres_restart 500  -vel_ksp_rtol 1.e-8 -vel_pc_type gamg -vel_ksp_converged_reason \
            -pres_ksp_type cg -pres_pc_type gamg -pres_ksp_rtol 1.e-4 -pres_ksp_converged_reason -ksp_atol 0.0"
end
    


t0 =0.0
dt = 0.002
tF = 15.0

Re = 60_000
D = 3
backend =  with_debug
rank_partition = (80,1,1)


sprob = StabilizedProblem(VMS(1))
timep = TimeParameters(t0=t0,dt=dt,tF=tF, t_endramp=2.0)

physicalp = PhysicalParameters(Re=Re)
solverp = SolverParameters(M=25, Number_Skip_Expansion=1e5, petsc_options = petsc_options_cstm(), matrix_freq_update=5)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["ph","friction"]])


meshp= MeshParameters(rank_partition,D,mesh_file)

simparams = SimulationParameters(timep,physicalp,solverp,exportp)


AirfoilCase = Airfoil(meshp,simparams,sprob)


SegregatedVMSSolver.solve(AirfoilCase,backend)

