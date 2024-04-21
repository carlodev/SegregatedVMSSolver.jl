module CommonsTestsSequential
using SegregatedVMSSolver
using Test
using PartitionedArrays

# include(joinpath("..","CommonsTests", "InitializeParamsTests.jl"))
# include(joinpath("..","CommonsTests","AddNewTagsTests.jl"))
# include(joinpath("..","CommonsTests","StabParamsTests.jl"))
# include(joinpath("..","CommonsTests","LinearUtilitiesTests.jl"))
# include(joinpath("..","CommonsTests","StabilizedEquationsTests.jl"))
include(joinpath("..","CommonsTests","MatrixCreationTests.jl"))


include(joinpath("..","CommonsTests", "ParametersDefTests.jl"))
include(joinpath("..","CommonsTests","AddNewTagsTests.jl"))
include(joinpath("..","CommonsTests","ModelCreationTests.jl"))
include(joinpath("..","CommonsTests","BoundaryConditionsTests.jl"))
include(joinpath("..","CommonsTests","SpaceConditionsTests.jl"))
include(joinpath("..","CommonsTests","EquationsTests.jl"))
include(joinpath("..","CommonsTests","RestartTests.jl"))
include(joinpath("..","CommonsTests","InitialConditionsTests.jl"))

function test_common_debug()
  with_debug() do distribute
  # ParametersDefTests.main(distribute)
  # AddNewTagsTests.main(distribute)
  # ModelCreationTests.main(distribute)
  # BoundaryConditionsTests.main(distribute)
  # SpaceConditionsTests.main(distribute)
  # EquationsTests.main(distribute)
  # RestartTests.main(distribute)
  # InitialConditionsTests.main(distribute)
  # MatrixCreationTests.main(distribute)

  # StabParamsTests.main(distribute)
  # LinearUtilitiesTests.main(distribute)
  # StabilizedEquationsTests.main(distribute)

  end
end

#with_debug
@testset "Commons Debug" begin
  test_common_debug()
end


end

