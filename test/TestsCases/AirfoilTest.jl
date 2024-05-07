using PartitionedArrays
using SegregatedVMSSolver
using Test
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions

function airfoil_test(backend)

t0 =0.0
dt = 1e-3
tF = 5e-3

Re = 10
D = 2
rank_partition = (2,2)
airfoil_mesh_file = joinpath(@__DIR__,"..", "..", "models", "DU89_2D_A1_M.msh")
airfoil_restart_file = joinpath(@__DIR__,"..", "..", "restarts", "BL_DU89_2D_A1_M.csv")



sprob = StabilizedProblem(VMS(1))
timep = TimeParameters(t0,dt,tF)

physicalp = PhysicalParameters(Re=Re)
solverp = SolverParameters(M=2)
exportp = ExportParameters(printinitial=false,printmodel=false,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])


meshp= MeshParameters(rank_partition,D,airfoil_mesh_file)
restartp = RestartParameters(airfoil_restart_file)

simparams = SimulationParameters(timep,physicalp,solverp,exportp,restartp)


mcase = Airfoil(meshp,simparams,sprob)


@test SegregatedVMSSolver.main(mcase,backend)
end

#airfoil_test(with_mpi)

#mpiexecjl -n 4 julia --project=. test/TestsCases/AirfoilTest.jl

