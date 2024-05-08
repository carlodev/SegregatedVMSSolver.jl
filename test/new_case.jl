
using PartitionedArrays
using Gridap
using GridapGmsh
using Revise
using Parameters
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using MPI

### test
t0 =0.0
dt = 0.1
tF = 20.0
t_endramp = 1.0

Re = 1000
D = 2
backend = with_mpi
rank_partition = (2,2)

# mesh_file = joinpath(@__DIR__, "..", "models", "DU89_2D_A1_M.msh")

sprob = StabilizedProblem(SUPG())
timep = TimeParameters(t0=t0,dt=dt,tF=tF)

physicalp = PhysicalParameters(Re=Re)
solverp = SolverParameters(petsc_options="")
exportp = ExportParameters(printinitial=false)

meshp= MeshParameters(rank_partition,D; N=50,L=0.5)

simparams = SimulationParameters(timep,physicalp,solverp,exportp)


mcase = LidDriven(meshp,simparams,sprob)

using SegregatedVMSSolver

SegregatedVMSSolver.main(mcase,backend)


#mpiexecjl --project=. -n 4 julia test/new_case.jl






# timep = TimeParameters(0.0,0.01,2.5)
# solverp = SolverParameters(petsc_options=" ", matrix_freq_update=1)

# physicalp = PhysicalParameters(Re=Re)

# exportp = ExportParameters()
# meshp= MeshParameters(rank_partition,D;N=50,L=0.5)
# simparams = SimulationParameters(timep,physicalp,solverp,exportp)

# mcase =TaylorGreen(meshp,simparams,sprob)
# using SegregatedVMSSolver

# SegregatedVMSSolver.main(mcase,with_debug)

