module ModelCreationTests


using Gridap
using Test
using Gridap
using GridapDistributed
using PartitionedArrays

using SegregatedVMSSolver
using SegregatedVMSSolver.CreateProblem
using SegregatedVMSSolver.ParametersDef



include(joinpath("..","..","case_test.jl")) 



function test_modelcreation(rank_partition, distribute, D)
    parts  = distribute(LinearIndices((prod(rank_partition),)))
    
    
    for (TestCase,mesh_file) in iterate_test_cases(D)
        simcase = create_simulation_test(TestCase, D; meshfile = mesh_file)
        model = create_model(parts, simcase)
    end

    return true
end



function main(distribute)
    for D in [2, 3]
        rank_partition = (D == 2) ? (2, 2) : (2, 2, 1)
        @test test_modelcreation(rank_partition, distribute, D)
    end
end

end
