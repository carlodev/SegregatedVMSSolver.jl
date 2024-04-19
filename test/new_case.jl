
using PartitionedArrays
using Gridap
using GridapGmsh
using Revise
using SegregatedVMSSolver


### test
t0 =0.0
dt = 0.1
tF = 1.0
Re = 1000
D = 2
backend = with_debug
rank_partition = (2,2)

mesh_file = joinpath(@__DIR__, "..", "models", "DU89_2D_A1_M.msh")

sprob = StabilizedProblem()
timep = TimeParameters(t0,dt,tF)
physicalp = PhysicalParameters(Re=Re)
solverp = SolverParameters(petsc_options=" ")
exportp = ExportParameters()
meshp= MeshParameters(rank_partition,D,mesh_file)


simparams = SimulationParameters(timep,physicalp,solverp,exportp)


mcase = Airfoil(meshp,simparams,sprob)


meshp= MeshParameters(rank_partition,D;N=10,L=0.5)
mcase =TaylorGreen(meshp,simparams,sprob)



using SegregatedVMSSolver
SegregatedVMSSolver.main(mcase,with_debug)

