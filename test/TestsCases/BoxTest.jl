using PartitionedArrays
using SegregatedVMSSolver

using SyntheticEddyMethod
using Test
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using SegregatedVMSSolver.CreateProblem


create_new_case(:Box)

function SegregatedVMSSolver.CreateProblem.create_boundary_conditions(simcase::Box, u_free,u_SEM, u_wall) 
        u_diri_tags=["inlet"]
        u_diri_values = [u_SEM]
        p_diri_tags=["outlet","limits"]
        p_diri_values = [0.0, 0.0]
        return u_diri_tags,u_diri_values,p_diri_tags,p_diri_values
end 



t0 =0.0
dt = 0.01

tF = 1.0

Re = 100
D = 3
##Turbulence Intensity
TI = 0.01

rank_partition = (2,2,1)

airfoil_mesh_file = joinpath(@__DIR__,"..", "..", "models", "Block_SEM.msh")


sprob = StabilizedProblem(VMS(1))
timep = TimeParameters(t0=t0,dt=dt,tF=tF)

physicalp = PhysicalParameters(Re=Re)

#Turbulence
Vboxinfo = VirtualBox((0,1.0), (0,1); σ=0.05)
Vboxinfo.N

# Vboxinfo.N=100 #reducing number of eddies


turbulencep= TurbulenceParameters(TI, Vboxinfo, physicalp)
turbulencep.Eddies[1]

solverp = SolverParameters(M=5,Number_Skip_Expansion=2)
# exportp = ExportParameters(printinitial=false,printmodel=false,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])
exportp = ExportParameters(printinitial=true,printmodel=true)

meshp= MeshParameters(rank_partition,D,airfoil_mesh_file)

simparams = SimulationParameters(timep,physicalp,turbulencep,solverp,exportp)


mcase = Box(meshp,simparams,sprob)

backend = with_mpi
SegregatedVMSSolver.solve(mcase,backend) 


# end

#turbulent_airfoil_test(with_mpi)
# turbulent_airfoil_test(with_debug)

#mpiexecjl -n 4 julia --project=. test/TestsCases/BoxTest.jl

