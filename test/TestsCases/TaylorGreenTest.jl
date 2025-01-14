using PartitionedArrays
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using MPI

function taylorgreen_test(backend)

t0 =0.0
dt = 0.1
tF = 0.5

Re = 1000
D = 2
rank_partition = (2,2)


sprob = StabilizedProblem(VMS(1))
timep = TimeParameters(t0=t0,dt=dt,tF=tF)

physicalp = PhysicalParameters(Re=Re)
solverp = SolverParameters()
exportp = ExportParameters(printinitial=false,printmodel=false)


meshp= MeshParameters(rank_partition,D;N=32,L=0.5)
simparams = SimulationParameters(timep,physicalp,solverp,exportp)



bc_tgv = Periodic(meshp,physicalp ) 

mcase = TaylorGreen(bc_tgv, meshp,simparams,sprob)


@test SegregatedVMSSolver.solve(mcase,backend)

end