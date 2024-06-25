using PartitionedArrays
using SegregatedVMSSolver
using SyntheticEddyMethod
using Test
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions

function turbulent_airfoil_test(backend)

t0 =0.0
dt = 1e-3
tF = 10e-3

Re = 10
D = 3
##Turbulence Intensity
TI = 0.001

rank_partition = (2,2,1)
airfoil_mesh_file = joinpath(@__DIR__,"..", "..", "models", "sd7003s_3D_simple.msh")


sprob = StabilizedProblem(VMS(1))
timep = TimeParameters(t0=t0,dt=dt,tF=tF, time_window=(2*dt,4*dt))

physicalp = PhysicalParameters(Re=Re)

#Turbulence
Vboxinfo = VirtualBox((-1,1), (0.0,0.20); σ=0.0125)

Vboxinfo.N=100 #reducing number of eddies

turbulencep= TurbulenceParameters(TI, Vboxinfo, physicalp)


solverp = SolverParameters(M=2,Number_Skip_Expansion=2)
exportp = ExportParameters(printinitial=false,printmodel=false,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])
meshp= MeshParameters(rank_partition,D,airfoil_mesh_file)


simparams = SimulationParameters(timep,physicalp,turbulencep,solverp,exportp)


mcase = Airfoil(meshp,simparams,sprob)


SegregatedVMSSolver.solve(mcase,backend) 

end

#turbulent_airfoil_test(with_mpi)
# turbulent_airfoil_test(with_debug)

#mpiexecjl -n 4 julia --project=. test/TestsCases/TurbulentAirfoilTest.jl

