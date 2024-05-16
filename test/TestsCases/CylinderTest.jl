using PartitionedArrays
using SegregatedVMSSolver
using Test
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions


function cylinder_test(backend)

t0 =0.0
dt = 0.1
t_endramp = 0.3
tF = 0.5

Re = 100
D = 2
rank_partition = (2,2)
cylinder_mesh_file = joinpath(@__DIR__,"..", "..", "models", "Cylinder_2D.msh")


sprob = StabilizedProblem(SUPG(1))
timep = TimeParameters(t0=t0,dt=dt,tF=tF,t_endramp=t_endramp)

physicalp = PhysicalParameters(Re=Re)
solverp = SolverParameters()
exportp = ExportParameters(printinitial=false,printmodel=false)


meshp= MeshParameters(rank_partition,D,cylinder_mesh_file)
simparams = SimulationParameters(timep,physicalp,solverp,exportp)


mcase = Cylinder(meshp,simparams,sprob)

@test SegregatedVMSSolver.solve(mcase,backend)

end