using SegregatedVMSSolver
using Parameters,PartitionedArrays
using Revise
using Test


@testset "Commons" begin include("CommonsTests/CommonsTests.jl") end


@testset "Case Test" begin include("case_test.jl") 
    run_case_test("TaylorGreen")
    run_case_test("LidDriven"; Reynolds = 10000)
    run_case_test("Cylinder"; meshfile="Cylinder_2D.msh", Reynolds = 1000)
end

@testset "Airfoil" begin include("airfoil_test.jl") 
    run_airfoil_test()
end

@testset "Utils" begin include("UtilsTests/UtilsTests.jl") end

