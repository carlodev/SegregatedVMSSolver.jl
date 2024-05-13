using PartitionedArrays
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions


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

log_dir = "Log"
mkdir(log_dir)
open(joinpath(log_dir,"PrintSim.txt"), "w") do io
end

mcase = TaylorGreen(meshp,simparams,sprob)


@test SegregatedVMSSolver.main(mcase,backend)

rm(log_dir, recursive=true)

end