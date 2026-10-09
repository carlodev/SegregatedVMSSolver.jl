using PartitionedArrays
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using MPI





t0 = 0.0
dt = 0.01
tF = 10*dt
vortex_diameter = 1.0
N = 16
Re = 1600
D = 2

backend = with_debug

rank_partition = (2,2)





# solver_options = petsc_options(; vel_ksp="gmres", vel_pc="gamg", pres_ksp="cg", pres_pc="gamg")
solver_options = " -vel_ksp_type gmres -vel_pc_type gamg -vel_ksp_rtol 1.e-10 -vel_ksp_converged_reason \
  -pres_ksp_type cg -pres_pc_type gamg -pres_ksp_rtol 1.e-6 -pres_ksp_converged_reason \
  -ksp_atol 0.0 -vec_type cuda -mat_type aijcusparse -log_view"
sprob = StabilizedProblem(VMS(3))
timep = TimeParameters(t0=t0, dt=dt, tF=tF)

physicalp = PhysicalParameters(Re=Re, c=vortex_diameter)
solverp = SolverParameters(matrix_freq_update=1, Number_Skip_Expansion=10e6, M=40,
petsc_options=solver_options)
exportp = ExportParameters(printinitial=true, printmodel=true, 
vtu_export = ["uh","ph","uh_analytic", "ph_analytic"], extra_export=["VelocityError","PressureError"])




meshp = MeshParameters(rank_partition, D; N=N, L=0.5 * vortex_diameter)
simparams = SimulationParameters(timep, physicalp, solverp, exportp)
bc_tgv = Periodic(meshp, physicalp)



mcase = TaylorGreen(bc_tgv, meshp, simparams, sprob)



# Create folder and file
# mkdir("Log")
# open("Log/PrintSim.txt", "w") do file
# end



SegregatedVMSSolver.solve(mcase, backend)

