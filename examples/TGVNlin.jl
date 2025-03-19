using PartitionedArrays
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions

using MPI





t0 = 0.0
dt = 0.01
tF = 2.5
vortex_diameter = 1.0
N = 16
Re = 1600
D = 2

backend = with_debug

rank_partition = (2,2)





solver_options =  "-snes_type newtonls -snes_linesearch_type basic -snes_linesearch_damping 1.0 -snes_rtol 1.0e-8 -snes_atol 0 -snes_monitor  -snes_max_it 20 \
-pc_type gamg -ksp_type gmres -ksp_converged_reason -ksp_max_it 50 -ksp_rtol 1e-8 -ksp_atol 0.0"

sprob = StabilizedProblem(VMS(1))
timep = TimeParameters(t0=t0, dt=dt, tF=tF)

physicalp = PhysicalParameters(Re=Re, c=vortex_diameter)
solverp = SolverParameters(matrix_freq_update=1, Number_Skip_Expansion=10e6, M=40,
petsc_options=solver_options, linear=true, θ=1.0)
exportp = ExportParameters(printinitial=false, printmodel=false, 
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

