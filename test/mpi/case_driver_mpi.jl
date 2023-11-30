using SegregatedVMSSolver
using MPI
using Test
using PartitionedArrays

include(joinpath("..","case_test.jl")) 
include(joinpath("..","airfoil_test.jl")) 


function testmpi()
    @testset "Case MPI" begin

            run_case_test("TaylorGreen", with_mpi)
            run_case_test("LidDriven",with_mpi; Reynolds = 10000)
            run_case_test("Cylinder",with_mpi; meshfile="Cylinder_2D.msh", Reynolds = 1000)
           
    end
    @testset "Airfoil MPI" begin

        run_airfoil_test(with_mpi)
end
end


testmpi()