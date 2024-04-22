
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.ModelCreation
using SegregatedVMSSolver.BoundaryConditions
using SegregatedVMSSolver.SpaceConditions
using SegregatedVMSSolver.SolveProblem


function main(simcase::SimulationCase,backend::Function)
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

        run_function(simcase,distribute)
    end

return true

end


function run_function(simcase::SimulationCase,distribute)
    params = Dict{Symbol,Any}()
    @sunpack rank_partition,order = simcase

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

    
    new_dict = Dict(:parts=>parts, :model=>model,
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
    
    solve_case(params,simcase)


end



