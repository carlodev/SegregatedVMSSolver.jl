




function liddriven_test()

t0 =0.0
dt = 0.1
tF = dt*5
t_endramp=0.3

Re = 1000
D = 2
rank_partition = (2,2)


sprob = StabilizedProblem(method=VMS(), coeff_method=ScalarFormulation())
timep = TimeParameters(t0=t0,dt=dt,tF=tF, t_endramp=t_endramp)

physicalp = PhysicalParameters(Re=Re)
solverp = SolverParameters()
exportp = ExportParameters(printinitial=true,printmodel=true)


meshp= MeshParameters(rank_partition,D;N=32,L=0.5)
simparams = SimulationParameters(timep,physicalp,solverp,exportp)


mcase = LidDriven(meshp,simparams,sprob)

log_dir = "Log"
mkpath(log_dir)
open(joinpath(log_dir,"PrintSim.txt"), "w") do io
end

@test typeof(mcase) <: SimulationCase

return mcase


end