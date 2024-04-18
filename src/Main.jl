function main(simcase::SimulationCase)
    #check(simcase)

    backend = get_field(simcase,:backend)
    
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