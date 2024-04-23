module RestartTests

using Test
using SegregatedVMSSolver
using Gridap
using GridapDistributed
using PartitionedArrays
using CSV
using DataFrames
using NearestNeighbors

using SegregatedVMSSolver.CreateProblem
using SegregatedVMSSolver.ParametersDef



include(joinpath("..","..","case_test.jl")) 


function test_restart(rank_partition, distribute, D)
    parts  = distribute(LinearIndices((prod(rank_partition),)))
    
    params = Dict{Symbol,Any}()


    t0 =0.0
    dt = 0.1
    tF = 1.0
    Re = 1000
    airfoil_mesh_file = joinpath(@__DIR__, "..","..","..", "models", "DU89_2D_A1_M.msh")
    airfoil_restart_file = joinpath(@__DIR__, "..", "..","..","restarts", "BL_DU89_2D_A1_M.csv")

    sprob = StabilizedProblem()
    timep = TimeParameters(t0,dt,tF)
    physicalp = PhysicalParameters(Re=Re)
    solverp = SolverParameters()
    exportp = ExportParameters(printmodel=false)
    restartp = RestartParameters(airfoil_restart_file)
    meshp= MeshParameters(rank_partition,D,airfoil_mesh_file)

    
    simparams = SimulationParameters(timep,physicalp,solverp,exportp,restartp)

    simcase = Airfoil(meshp,simparams,sprob)
    @sunpack order = simcase

        model = create_model(parts, simcase)
    
        boundary_conditions = create_boundary_conditions(simcase) 
    
        V, U, P, Q, Y, X = creation_fe_spaces(simcase, model, boundary_conditions)
       

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
        
        @sunpack restartfile,D = simcase
        restart_df = DataFrame(CSV.File(restartfile))

        tree = SegregatedVMSSolver.CreateProblem.create_search_tree(restart_df)
        @test typeof(tree)<:BruteTree

        uh = SegregatedVMSSolver.CreateProblem.restart_uh_field(D,tree,restart_df)
        @test typeof(uh)<:Function

        ph = SegregatedVMSSolver.CreateProblem.restart_ph_field(tree,restart_df)
        @test typeof(ph)<:Function


        return true

end



function main(distribute)
        D = 2
        rank_partition = (2, 2)
        @test test_restart(rank_partition, distribute, D)
end





end