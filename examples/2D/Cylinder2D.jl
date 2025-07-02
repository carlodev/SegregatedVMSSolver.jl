using PartitionedArrays
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using MPI





t0 = 0.0
dt = 0.1
tF = 10.0
t_endramp=2.0

Re = 1000
D = 2

backend = with_debug

cylinder_mesh_file = joinpath(@__DIR__,"..","..",  "models", "Cylinder_2D.msh")


rank_partition = (2,2)



solver_options = " -vel_ksp_type gmres -vel_pc_type gamg -vel_ksp_rtol 1.e-4 -vel_ksp_converged_reason \
    -pres_ksp_type cg -pres_pc_type gamg  -pres_ksp_rtol 1.e-2 -pres_ksp_converged_reason \
    -ksp_atol 0.0"

sprob = StabilizedProblem(VMS(2))
timep = TimeParameters(t0=t0, dt=dt, tF=tF, t_endramp=t_endramp)

physicalp = PhysicalParameters(Re=Re, c=1.0)
solverp = SolverParameters(matrix_freq_update=1, Number_Skip_Expansion=10e6, M=4,
    petsc_options=solver_options)
exportp = ExportParameters(printinitial=true, printmodel=true)


meshp= MeshParameters(rank_partition,D,cylinder_mesh_file)


simparams = SimulationParameters(timep, physicalp, solverp, exportp)

mcase = Cylinder(meshp,simparams,sprob)



# Create folder and file
mkdir("Log")
open("Log/PrintSim.txt", "w") do file
end


SegregatedVMSSolver.solve(mcase, backend)









