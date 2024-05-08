module BoundaryConditionsTests

using Gridap
using Test
using Gridap
using GridapDistributed
using PartitionedArrays

using SegregatedVMSSolver
using SegregatedVMSSolver.CreateProblem
using SegregatedVMSSolver.ParametersDef



include(joinpath("..","..","case_test.jl")) 



function test_boundaryconditions(distribute, D)
    
    for (TestCase,mesh_file) in iterate_test_cases(D)
        simcase = create_simulation_test(TestCase, D; meshfile = mesh_file)
        boundary_conditions = create_boundary_conditions(simcase) 
        @test typeof(boundary_conditions) <:Tuple
    end
end


function main(distribute)
    for D in [2, 3]
        test_boundaryconditions(distribute, D)
    end
end

end
