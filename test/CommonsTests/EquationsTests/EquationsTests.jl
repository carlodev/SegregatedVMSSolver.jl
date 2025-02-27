module EquationsTests

using Test
using SegregatedVMSSolver
using Gridap
using Gridap.Fields
using GridapDistributed
using PartitionedArrays

using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.CreateProblem
using SegregatedVMSSolver.Equations



include("StabilizedEquationsTests.jl")
include("StabilizationParametersTests.jl")
include("StabilizationOperationsTests.jl")


function main(uh, params, mcase)
    @testset "Equations Tests" begin
        D = mcase.meshp.D
        @test StabilizedEquationsTests.main(uh, params, mcase)
        @test StabilizationParametersTests.main(params[:Ω], D)
        @test StabilizationOperationTests.main(mcase, params, uh, params[:Ω], D)

    end
end


end #end module