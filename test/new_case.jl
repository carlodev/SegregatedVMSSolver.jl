
using PartitionedArrays
using Gridap
using GridapGmsh

### test
t0 =0.0
dt = 0.1
tF = 1.0
Re = 1000
D = 2
backend = with_debug
rank_partition = (2,2)

mesh_file = joinpath(@__DIR__, "..", "models", "DU89_2D_A1_M.msh")

sprob = StabilizedProblem()
timep = TimeParameters(t0,dt,tF)
physicalp = PhysicalParameters(Re=Re,D=D)
solverp = SolverParameters()
exportp = ExportParameters()
clusterp= ClusterParameters(rank_partition,backend)
simparams = SimulationParameters(timep,physicalp,solverp,exportp,clusterp)

airfoil_case = Airfoil(mesh_file,simparams,sprob)

printstructure(airfoil_case)

new_main(airfoil_case)



function run_function(simcase,distribute)
    params = Dict{Symbol,Any}()
    rank_partition,order = get_field(simcase,[:rank_partition,:order])

    parts  = distribute(LinearIndices((prod(rank_partition),)))

    model = create_model(parts, simcase)
    @info "model read completed"

    boundary_conditions = create_boundary_conditions(simcase) 
    @info "boundary conditions created"

    V, U, P, Q, Y, X = creation_fe_spaces(simcase, model, boundary_conditions)
    @info "FE Spaces Created"

    trials = [U, P]
    tests = [V, Q]
    
    degree = 4*order
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    
    new_dict = Dict(:parts=>parts,
    :U => U,
    :P => P,
    :X => X,
    :Y => Y,
    :Ω => Ω,
    :dΩ => dΩ,
    :degree => degree,
    :trials => trials, 
    :tests => tests)
    merge!(params, new_dict)

  
    # solve_case(params,simcase)
end


