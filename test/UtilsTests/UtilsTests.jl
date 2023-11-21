module Utilstests

using Test

@testset "Utils" begin
  include("ReadAirfoilResultsTests.jl")

  include("FindWallDistanceTests.jl")

end


end