using PartitionedArrays

using SegregatedVMSSolver

using SyntheticEddyMethod
using Test
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions



function turbulent_airfoil_test(backend, create_sem_boundary)

t0 =0.0

dt = 0.001

tF = dt * 2

Re = 1000
D = 2

##Turbulence Intensity
TI = 0.1

rank_partition = (2,2)
airfoil_mesh_file = joinpath(@__DIR__,"..", "..", "models", "DU89_2D_A1_M.msh")


sprob = StabilizedProblem(VMS(1))
timep = TimeParameters(t0=t0,dt=dt,tF=tF)

physicalp = PhysicalParameters(Re=Re)
initialp = InitialParameters(u0=[1.0,0.0])

#Turbulence
Vboxinfo = VirtualBox((-1,1), (0.0,0.2); σ=0.0125)
@test typeof(Vboxinfo.N) == Int64


turbulencep= TurbulenceParameters(TI, Vboxinfo, physicalp, create_sem_boundary)


solverp = SolverParameters(M=2,Number_Skip_Expansion=2)
# exportp = ExportParameters(printinitial=true,printmodel=true)
exportp = ExportParameters(printinitial=false,printmodel=false)


meshp= MeshParameters(rank_partition,D,airfoil_mesh_file)

simparams = SimulationParameters(timep,physicalp,turbulencep,solverp,exportp,initialp)


mcase = Airfoil(meshp,simparams,sprob)

@test SegregatedVMSSolver.solve(mcase,backend) 

end


#mpiexecjl -n 4 julia --project=. test/TestsCases/TurbulentAirfoilTest.jl

