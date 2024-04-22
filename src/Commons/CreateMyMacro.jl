
using SegregatedVMSSolver
using Parameters
using SegregatedVMSSolver.ParametersDef


timep = TimeParameters(0.0,0.5,1.0)
physicalp = PhysicalParameters(Re=1000)
solverp = SolverParameters(petsc_options="")
exportp = ExportParameters(printinitial=false)
meshp= MeshParameters((2,2),2; N=50,L=0.5)

simparams = SimulationParameters(timep,physicalp,solverp,exportp)
mcase = LidDriven(meshp,simparams,sprob)


@sunpack D,rank_partition,order,L,θ,log_dir,order,method= mcase