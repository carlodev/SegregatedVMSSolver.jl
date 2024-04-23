module InitialConditionsTests

using Test
using SegregatedVMSSolver
using Gridap
using GridapDistributed
using PartitionedArrays

using SegregatedVMSSolver.CreateProblem
using SegregatedVMSSolver.ParametersDef



include(joinpath("..","..","case_test.jl")) 



function test_initialconditions(rank_partition, distribute, D)
    parts  = distribute(LinearIndices((prod(rank_partition),)))
    
    params = Dict{Symbol,Any}()

    TestCase,mesh_file = first(iterate_test_cases(D))

    simcase = create_simulation_test(TestCase, D; meshfile = mesh_file)
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
        
        uh0,ph0 = create_initial_conditions(simcase,params)

        return true

end



function main(distribute)
        D = 2
        rank_partition = (2, 2)

        @test test_initialconditions(rank_partition, distribute, D)
end





end