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
D = 3

backend = with_debug

rank_partition = (1,1,1)




sol_options= petsc_options(; vel_ksp="gmres", vel_pc="gamg", pres_ksp = "cg", pres_pc = "asm")


sprob = StabilizedProblem(VMS(2))
timep = TimeParameters(t0=t0, dt=dt, tF=tF)

physicalp = PhysicalParameters(Re=Re, c=vortex_diameter)
solverp = SolverParameters(matrix_freq_update=1, Number_Skip_Expansion=10e6, M=40,petsc_options=sol_options)
exportp = ExportParameters(printinitial=true, printmodel=true, 
vtu_export = ["uh","ph"], extra_export=["KineticEnergy", "Enstrophy"])

meshp = MeshParameters(rank_partition, D; N=N, L=pi * vortex_diameter)


simparams = SimulationParameters(timep, physicalp, solverp, exportp)

bc_tgv = Periodic(meshp, physicalp)



mcase = TaylorGreen(bc_tgv, meshp, simparams, sprob)



# # Create folder and file
# mkdir("Log")
# open("Log/PrintSim.txt", "w") do file
# end



SegregatedVMSSolver.solve(mcase, backend)









