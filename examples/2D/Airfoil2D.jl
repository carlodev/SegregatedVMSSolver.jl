using PartitionedArrays
using SegregatedVMSSolver
using Test
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions


t0 =0.0
dt = 0.01
tF = 5.0

Re = 1000
D = 2
rank_partition = (2,2)
airfoil_mesh_file = joinpath(@__DIR__,"..", "..", "models", "DU89_2D_A1_M.msh")
airfoil_restart_file = joinpath(@__DIR__,"..", "..","restarts", "BL_DU89_2D_A1_M.csv")

backend = with_debug


sprob = StabilizedProblem(VMS(1))
timep = TimeParameters(t0=t0,dt=dt,tF=tF)


solver_options = " -vel_ksp_type gmres -vel_pc_type gamg -vel_ksp_rtol 1.e-4 -vel_ksp_converged_reason \
    -pres_ksp_type cg -pres_pc_type lu  -pres_ksp_rtol 1.e-2 -pres_ksp_converged_reason \
    -ksp_atol 0.0"

physicalp = PhysicalParameters(Re=Re)
solverp = SolverParameters(matrix_freq_update=10,M=10,Number_Skip_Expansion=1e5,petsc_options=solver_options)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])


meshp= MeshParameters(rank_partition,D,airfoil_mesh_file)
intialp = InitialParameters(airfoil_restart_file)

simparams = SimulationParameters(timep,physicalp,solverp,exportp,intialp)


mcase = Airfoil(meshp,simparams,sprob)


# Create folder and file
mkdir("Log")
open("Log/PrintSim.txt", "w") do file
end


SegregatedVMSSolver.solve(mcase, backend)









