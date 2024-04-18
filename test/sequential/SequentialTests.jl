module SequentialTests
using Test
using PartitionedArrays

@testset "Commons Sequential" begin include("CommonsTestsSequential.jl") end

@testset "Case Test" begin include(joinpath("..","case_test.jl")) 
    # run_case_test("TaylorGreen", with_debug)
    run_case_test("LidDriven",with_debug; Reynolds = 10000)
    # run_case_test("Cylinder",with_debug; meshfile="Cylinder_2D.msh", Reynolds = 1000)
end

# @testset "Airfoil" begin include(joinpath("..","airfoil_test.jl"))
#     run_airfoil_test(with_debug)
# end


end