module CommonsTestsSequential
using SegregatedVMSSolver
using Test
using PartitionedArrays
include(joinpath("..","CommonsTests", "InitializeParamsTests.jl"))
include(joinpath("..","CommonsTests","AddNewTagsTests.jl"))
include(joinpath("..","CommonsTests","StabParamsTests.jl"))
include(joinpath("..","CommonsTests","LinearUtilitiesTests.jl"))
include(joinpath("..","CommonsTests","StabilizedEquationsTests.jl"))
include(joinpath("..","CommonsTests","MatrixCreationTests.jl"))
include(joinpath("..","CommonsTests","RestartTests.jl"))

function test_common_debug()
  with_debug() do distribute
  InitializeParamsTests.main(distribute)
  AddNewTagsTests.main(distribute)
  StabParamsTests.main(distribute)
  LinearUtilitiesTests.main(distribute)
  StabilizedEquationsTests.main(distribute)
  MatrixCreationTests.main(distribute)
  RestartTests.main(distribute)
  end
end

#with_debug
@testset "Commons Debug" begin
  test_common_debug()
end

end