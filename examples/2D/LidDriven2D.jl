using PartitionedArrays
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using MPI





t0 = 0.0
dt = 0.1
tF = 10.0
t_endramp=2.0

N = 16
Re = 1000
D = 2

backend = with_debug

rank_partition = (1,1)



solver_options = petsc_options(; vel_ksp="gmres", vel_pc="gamg", pres_ksp="cg", pres_pc="gamg")

sprob = StabilizedProblem(VMS(2))
timep = TimeParameters(t0=t0, dt=dt, tF=tF, t_endramp=t_endramp)

physicalp = PhysicalParameters(Re=Re, c=1.0)
solverp = SolverParameters(matrix_freq_update=1, Number_Skip_Expansion=10e6, M=40,
    petsc_options=solver_options)
exportp = ExportParameters(printinitial=true, printmodel=true)


meshp = MeshParameters(rank_partition, D; N=N, L=0.5)


simparams = SimulationParameters(timep, physicalp, solverp, exportp)

mcase = LidDriven(meshp,simparams,sprob)




# Create folder and file
mkdir("Log")
open("Log/PrintSim.txt", "w") do file
end


SegregatedVMSSolver.solve(mcase, backend)









