module SpaceConditionsTests

using Gridap
using Test
using Gridap
using GridapDistributed
using PartitionedArrays

using SegregatedVMSSolver
using SegregatedVMSSolver.ModelCreation
using SegregatedVMSSolver.BoundaryConditions
using SegregatedVMSSolver.SpaceConditions



include(joinpath("..","case_test.jl")) 


function test_spaceconditions(rank_partition, distribute, D)
    parts  = distribute(LinearIndices((prod(rank_partition),)))
    

    for (TestCase,mesh_file) in iterate_test_cases(D)
        simcase = create_simulation_test(TestCase, D; meshfile = mesh_file)

        model = create_model(parts, simcase)
    
        boundary_conditions = create_boundary_conditions(simcase) 
    
        fespaces = creation_fe_spaces(simcase, model, boundary_conditions)
        @test typeof(fespaces)<:Tuple
        @test length(fespaces) == 6
    end

end



function main(distribute)
    for D in [2, 3]
        rank_partition = (D == 2) ? (2, 2) : (2, 2, 1)
        test_spaceconditions(rank_partition, distribute, D)
    end
end

end
