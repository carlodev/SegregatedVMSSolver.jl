using MPI
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.CreateProblem
using SegregatedVMSSolver.SolveProblem
using PartitionedArrays, Gridap, GridapDistributed

function solve(simcase::SimulationCase,backend::Function)
    #check(simcase)

    backend() do distribute
        if backend == with_mpi
            comm = MPI.COMM_WORLD
            #To avoid multiple printing of the same line in parallel
            if MPI.Comm_rank(comm) != 0
              redirect_stderr(devnull)
              redirect_stdout(devnull)
            end
           
        end

        printstructure(simcase)
        params = setup_case(simcase,distribute)
        solve_case(params,simcase)
    end

return true

end


function setup_case(simcase::SimulationCase,distribute)
    params = Dict{Symbol,Any}()
    @sunpack rank_partition,order = simcase

    parts  = distribute(LinearIndices((prod(rank_partition),)))

    model = create_model(parts, simcase)
    @info "model read completed"

    boundary_conditions = create_boundary_conditions(simcase) 
    @info "boundary conditions created"

    V, U, P, Q = creation_fe_spaces(simcase, model, boundary_conditions)

    @info "FE Spaces Created"

    trials = [U, P]
    tests = [V, Q]
    
    degree = 2*order
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    
    new_dict = Dict(:parts=>parts, :model=>model,
    :U=>U,:V=>V,:P=>P,:Q=>Q,
    :Ω => Ω,
    :dΩ => dΩ,
    :degree => degree,
    :trials => trials, 
    :tests => tests)
    merge!(params, new_dict)
    
    return params


end



