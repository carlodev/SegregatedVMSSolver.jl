module ParametersDefTests

using Test
using SegregatedVMSSolver
using Parameters
using PartitionedArrays


include(joinpath("..","case_test.jl")) 
 
function main(distribute)
    D = 2
    
    for (TestCase,mesh_file) in iterate_test_cases(D)
        tcase = create_simulation_test(TestCase, D; meshfile = mesh_file)
        @test typeof(tcase) <: SimulationCase
    end
            
end


end #end module